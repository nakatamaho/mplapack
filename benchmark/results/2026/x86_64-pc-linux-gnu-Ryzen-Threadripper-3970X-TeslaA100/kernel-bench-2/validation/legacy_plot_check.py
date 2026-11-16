"""Compare old/new legacy-default rendering of an actual historical CSV."""
import csv,hashlib,pathlib,subprocess,sys,tempfile,types
here=pathlib.Path(__file__).resolve()
sys.path.insert(0,str(here.parents[5]))
import plot_sweep
base=types.ModuleType('plot_sweep_before_item4')
code=subprocess.check_output(['git','show','fc85906fc:benchmark/plot_sweep.py'],text=True,cwd=here.parents[6])
exec(compile(code,'plot_sweep_before_item4.py','exec'),base.__dict__)
p=here.parents[2]/'rgemm/mpfr512_sweep10k.csv'
with p.open() as f:rows=[(int(r['n']),float(r['gflops'])) for r in csv.DictReader(f)]
with tempfile.TemporaryDirectory() as tmp:
    d=pathlib.Path(tmp)
    for dark in (False,True):
        style='dark' if dark else 'light'
        a=d/f'before-{style}.png';b=d/f'after-{style}.png'
        base.draw(rows,a,dark=dark);plot_sweep.draw(rows,b,dark=dark)
        same=a.read_bytes()==b.read_bytes()
        print(f'style={style} byte_identical={same} sha256={hashlib.sha256(b.read_bytes()).hexdigest()}')
        assert same
print('PASS: actual historical MPFR512 CSV, both styles, legacy defaults unchanged')
