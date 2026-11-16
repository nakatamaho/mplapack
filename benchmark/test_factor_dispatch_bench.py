#!/usr/bin/env python3
"""CLI/parser regressions, without CUDA or numerical timing claims."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class DispatchBenchTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.output = self.directory / 'samples.jsonl'
        self.runner = Path(__file__).with_name('run_factor_dispatch_bench.py')
        self.script = '''#!/usr/bin/env python3
import os, sys
from pathlib import Path
routine = next((x[2:] for x in sys.argv if x in ['--getri', '--potrf']), 'getrf')
name = ('C' if '--complex' in sys.argv else 'R') + routine
n = sys.argv[sys.argv.index('--n')+1]
print('routine,n,nrhs,call,temperature,seconds')
for i, seconds in enumerate([10, 2, 4]):
    print(f'{name},{n},-1,{i},{"cold" if i==0 else "warm"},{seconds}')
    if 'cuda_total' in sys.argv[0] and not Path(sys.argv[0]).with_name('missing').exists():
        print(f'{name} MPLAPACK_DD CUDA n={n} upload=0 inverse=0 download=0', file=sys.stderr)
'''
        for suffix in ['opt', 'cuda_total']:
            path = self.directory / ('Rfactor.dd_' + suffix)
            path.write_text(self.script)
            path.chmod(0o755)

    def run_cli(self, *extra):
        return subprocess.run([sys.executable, str(self.runner), '--bin-dir', str(self.directory),
                               '--backends', 'dd', '--sizes', '3', '--routines', 'getri',
                               '--kinds', 'real', '--rounds', '2', '--samples',
                               '--output', str(self.output), *extra], text=True, capture_output=True)

    def test_samples_and_resume(self):
        result = self.run_cli()
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = [json.loads(line) for line in self.output.read_text().splitlines()]
        self.assertEqual(len(rows), 6)
        self.assertEqual([r['mode'] for r in rows], ['cpu', 'default', 'forced', 'forced', 'default', 'cpu'])
        self.assertTrue(all(r['cold_seconds'] == 10 and r['seconds'] == 3 for r in rows))
        self.assertEqual(rows[2]['entry_trace_count'], 3)
        self.assertIn('cpu_cold,cpu_warm', result.stdout)
        self.assertEqual(self.run_cli().returncode, 0)
        self.assertEqual(len(self.output.read_text().splitlines()), 6)

    def test_missing_forced_trace(self):
        (self.directory / 'missing').touch()
        result = self.run_cli('--modes', 'forced')
        self.assertNotEqual(result.returncode, 0)
        row = json.loads(self.output.read_text())
        self.assertEqual(row['returncode'], 77)
        self.assertIn('forced entry rejected', result.stderr)

    def test_mpfr_inverse_not_claimed_as_gpu(self):
        result = self.run_cli('--backends', 'mpfr')
        self.assertEqual(result.returncode, 2)
        self.assertIn('MPFR inverse GPU entries do not exist', result.stderr)

    def test_samples_require_warm_call(self):
        self.assertEqual(self.run_cli('--repeat', '1').returncode, 2)


if __name__ == '__main__':
    unittest.main()
