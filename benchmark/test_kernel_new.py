#!/usr/bin/env python3
"""Host-only grid/budget/CLI smoke tests; no CUDA queries or GPU execution."""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import csv
from fractions import Fraction
from unittest.mock import patch, Mock
from types import SimpleNamespace
import kernel_new as kn

class Policy(unittest.TestCase):
    def test_vectors(self):
        v = kn.grid('Rdot')
        self.assertEqual(len(v), 57)
        self.assertEqual(v[0], 1024)
        self.assertEqual(v[-1], 1 << 24)
        self.assertTrue(all(1 << e in v for e in range(10,25)))
    def test_square_grid(self):
        for dense in (True,False):
            v = kn.grid('Rtrsv',dense)
            self.assertTrue(all(n % 7 == 0 or n % 31 == 0 or n & (n-1) == 0 for n in v))
            self.assertTrue(all(n in v for n in (1024,2048,4096)))
    def test_gpu_process_identity(self):
        with tempfile.TemporaryDirectory() as d:
            result=SimpleNamespace(returncode=0, stdout='42, Rkernel_bench.dd_cuda_kernel\n', stderr='')
            with patch('kernel_new.subprocess.run', return_value=result):
                self.assertTrue(kn.snapshot(Path(d)/'log','before',True)[1])
                self.assertFalse(kn.snapshot(Path(d)/'log','during',True,42)[1])
                self.assertTrue(kn.snapshot(Path(d)/'log','during',True,99)[1])
    def test_bundle(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d);(out/'raw').mkdir()
            row=dict.fromkeys(kn.FIELDS,'');row.update(type='dd',kernel='Rdot',method='cpu',n=1024,phase='pilot',exit=0)
            kn.append_csv(out/'kernels.csv',row)
            p=out/'raw/dd_Rdot_cpu_1024_pilot.log';p.write_text('unit test log\nexit=0\n')
            p.with_suffix('.log.exit').write_text('0\n')
            kn.bundle(SimpleNamespace(out=d))
            target=out/'raw/dd_Rdot_cpu.log';before=target.read_bytes()
            kn.bundle(SimpleNamespace(out=d))
            self.assertEqual(before,target.read_bytes())
            self.assertIn(b'exit receipt: 0',before)
            self.assertFalse(p.exists())
    def test_inverse_extension_stops(self):
        with tempfile.TemporaryDirectory() as d:
            r=kn.Runner.__new__(kn.Runner)
            r.a=SimpleNamespace(types='mpfr512',methods='gpu',kernels='Rtrtri',dense=False)
            r.out=Path(d);r.limit=1000;r.spent=0
            r.rows=[dict(type='mpfr512',kernel='Rtrtri',method='gpu',phase='sweep',n='512',seconds='1',elapsed='2')]
            r.measure=Mock(return_value=dict(n='1022',seconds='11',exit=0))
            r.extend_mpfr();self.assertEqual(r.measure.call_count,1)
            r.spent=999;r.measure.reset_mock();r.extend_mpfr()
            r.measure.assert_not_called()
            r.a.phase='extend-mpfr';r.a.kernels='RCnrm2';r.limit=65;r.spent=0
            r.rows=[dict(type='mpfr512',kernel='RCnrm2',method='gpu',phase='sweep',n=str(1<<23),seconds='8',elapsed='9')]
            r.measure=Mock(return_value=dict(n=str(5*(1<<21)),seconds='11',exit=0))
            r.extend_mpfr()
            r.measure.assert_called_once_with('mpfr512','RCnrm2','gpu',5*(1<<21),'sweep')
    def test_inventory(self):
        jobs = list(kn.series(['dd','td','qd','mpfr512'],['cpu','gpu']))
        self.assertEqual(len(jobs),180)
        self.assertEqual(len({(t,k) for t,k,m in jobs}),76)
        self.assertTrue(all(t == 'mpfr512' for t,k,m in jobs if k in kn.NAMES[17:]))


def smoke(bindir, typ):
    b = 'mpfr' if typ == 'mpfr512' else typ
    program = str(Path(bindir).resolve() / f'Rkernel_bench.{b}_opt')
    env = dict(os.environ, OMP_NUM_THREADS='32', OMP_DYNAMIC='FALSE', MPFRXX_DEFAULT_PRECISION_BITS='512', CUDA_VISIBLE_DEVICES='')
    kernels = kn.NAMES if b == 'mpfr' else kn.NAMES[:17]
    with tempfile.TemporaryDirectory() as d:
        for k in kernels:
            previous = None
            for mode in ('fast','default') if k in kn.SUMS or k in kn.NAMES[17:] else ('fast',):
                cmd = [program,k,mode,'7','--reps','1','--warmup','0','--dump',d+'/result']
                p = subprocess.run(cmd,env=env,capture_output=True,text=True)
                print('$',' '.join(cmd),'exit=',p.returncode, p.stdout,p.stderr,flush=True)
                assert p.returncode == 0 and Path(d+'/result').stat().st_size
                data = Path(d+'/result').read_bytes()
                if k in kn.NAMES[17:] and previous is not None:
                    assert previous == data, 'inverse must not opt its leaf reductions into FAST'
                previous = data
        for args in (['Rdot','invalid','7'],['Rdot','fast','0'],['Rdot','fast','1024','--reps','0']):
            p = subprocess.run([program]+args,env=env,capture_output=True,text=True)
            print('reject',args,'exit=',p.returncode,flush=True)
            assert p.returncode != 0

def exact_hex(text):
    # Sum BF limbs exactly; MPFR's %Ra has an arbitrary-length significand.
    total=Fraction(0)
    for limb in text.split(','):
        sign=-1 if limb.startswith('-') else 1
        h,e=limb.lstrip('+-').lower().split('p')
        whole,_,part=h.removeprefix('0x').partition('.')
        total += sign*Fraction(int(whole+part,16),16**len(part))*Fraction(2)**int(e)
    return total


def gpu_check(bindir, mpfr_bindir, out):
    if not kn.GATE.is_file():
        raise RuntimeError('GPU gate absent; refusing device checks')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    env=dict(os.environ, OMP_NUM_THREADS='32', OMP_DYNAMIC='FALSE',MPFRXX_DEFAULT_PRECISION_BITS='512')
    results=[]
    with tempfile.TemporaryDirectory() as d:
        for typ in ('dd','td','qd','mpfr512'):
            b='mpfr' if typ=='mpfr512' else typ
            base=Path(mpfr_bindir if b=='mpfr' else bindir).resolve()
            for k in kn.NAMES if b=='mpfr' else kn.NAMES[:17]:
                for n in (64,200):
                    dumps=[]
                    for method in ('cpu','gpu'):
                        suffix=('cuda_resident' if b=='mpfr' else 'cuda_kernel') if method=='gpu' else 'opt'
                        fn=Path(d)/method
                        cmd=[str(base/f'Rkernel_bench.{b}_{suffix}'),k,'fast',str(n),'--reps','3','--warmup','1','--dump',str(fn)]
                        p=subprocess.run(cmd,env=dict(env,OMP_NUM_THREADS='1' if method=='gpu' else '32'),capture_output=True,text=True)
                        log=out/f'{typ}_{k}_{n}_{method}.log'
                        log.write_text('$ '+' '.join(cmd)+'\n'+p.stdout+p.stderr+f'\nexit={p.returncode}\n')
                        log.with_suffix('.log.exit').write_text(str(p.returncode)+'\n')
                        if p.returncode:
                            raise RuntimeError(f'{cmd} exit={p.returncode}: {p.stderr}')
                        dumps.append(fn.read_text())
                    same=dumps[0]==dumps[1]
                    status='identical' if same else 'DIFFERENT'
                    max_error=Fraction(0)
                    if not same and k in kn.SUMS:
                        aa=dumps[0].strip().split(';');bb=dumps[1].strip().split(';')
                        u=Fraction(2)**(-{'dd':104,'td':157,'qd':209,'mpfr512':512}[typ])
                        passed=True
                        for a,bv in zip(aa,bb):
                            ref=exact_hex(a);err=abs(ref-exact_hex(bv));max_error=max(max_error,err)
                            bound=(16+8*(2*n-1).bit_length())*4*u*abs(ref) if 'nrm2' in k else 32*n*n*4*u
                            passed &= err<=bound
                        status='bounded-fast' if passed and len(aa)==len(bb) else 'DIFFERENT'
                    results.append([typ,k,n,status,float(max_error)])
                    print(*results[-1],flush=True)
    with (out/'checks.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['type','kernel','n','status','max_absolute_error']);w.writerows(results)
    if any(r[3]=='DIFFERENT' for r in results):
        raise RuntimeError('benchmark operand mismatch; see checks.csv')


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--bindir')
    p.add_argument('--type',default='dd')
    p.add_argument('--gpu-check',action='store_true')
    p.add_argument('--mpfr-bindir')
    p.add_argument('--out')
    a=p.parse_args()
    if a.gpu_check:
        gpu_check(a.bindir,a.mpfr_bindir,a.out)
    elif a.bindir:
        smoke(a.bindir,a.type)
    else:
        unittest.main(argv=[__file__])
