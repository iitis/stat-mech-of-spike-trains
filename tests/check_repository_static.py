"""Check repository wiring without MATLAB or third-party Python packages.

This validates paths and experiment configuration, not MATLAB execution.
Run: python tests/check_repository_static.py
"""
from pathlib import Path
import ast
import re
import shutil
import subprocess


def main():
    root = Path(__file__).resolve().parents[1]
    scripts = root / 'scripts/current'
    for folder in ['src', 'tests', 'scripts']:
        for path in (root / folder).rglob('*'):
            if path.suffix in {'.m', '.py', '.sh', '.sbatch'}:
                path.read_bytes().decode('utf-8')
    for doc in ['README.md', 'README_REPRODUCIBILITY.md']:
        text = (root / doc).read_text(encoding='utf-8')
        for reference in re.findall(r"run\('([^']+)'\)", text):
            assert (root / reference).is_file(), f'{doc}: missing {reference}'
    for path in (root / 'tests').glob('*.py'):
        ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
    ast.parse((scripts / 'plot_phase_multi.py').read_text(encoding='utf-8'))
    experiments = {
        'main': [('0_2', (0.0, 2.0, 0.2)), ('2p2_4', (2.2, 4.0, 0.2)), ('4p2_6', (4.2, 6.0, 0.2))],
        'convergence': [('0_2', (0.0, 2.0, 0.5)), ('2p5_4', (2.5, 4.0, 0.5)), ('4p5_6', (4.5, 6.0, 0.5))],
    }
    for kind, batches in experiments.items():
        for batch, beta in batches:
            path = scripts / f'run_{kind}_beta_{batch}.m'
            text = re.sub(r'%[^\n]*', '', path.read_text(encoding='utf-8'))
            args = text.split('run_phase_diagram(', 1)[1].split(');', 1)[0].replace('...', '')
            args = [a.strip() for a in args.split(',')]
            assert len(args) == 17, f'{path.name}: argument count'
            assert tuple(map(float,args[10:13])) == beta, path.name
            assert tuple(map(float,args[1:4])) == (100, 0.3, 5), path.name
            assert tuple(map(float,args[7:10])) == (500, 500, 3), path.name
            assert float(args[13]) == (100000 if kind == 'main' else 1000000), path.name
            assert float(args[14]) == (10000 if kind == 'main' else 100000), path.name
            assert tuple(map(float,args[4:7])) == (0, 4.5, 0.15 if kind == 'main' else 0.9), path.name
            assert args[15] == "'heaviside3'" and float(args[16]) == 1, path.name
            if kind == 'main':
                wrapper = (scripts / f'run_main_beta_{batch}.sbatch').read_text(encoding='utf-8')
                assert f'run_main_beta_{batch}.m' in wrapper and 'run_convergence_' not in wrapper
    for path in scripts.glob('*.sbatch'):
        for reference in re.findall(r"run\('([^']+)'\)", path.read_text(encoding='utf-8')):
            assert (root / reference).is_file(), f'{path.name}: missing {reference}'
    test_runner = (root / 'tests/run_all_tests.m').read_text(encoding='utf-8')
    for name in re.findall(r'^(test_\w+);', test_runner, re.MULTILINE):
        assert (root / f'tests/{name}.m').is_file(), name
    merge = (scripts / 'merge_missing_exprnd_points.m').read_text(encoding='utf-8')
    assert 'save(base_file' not in merge
    assert 'seed_grid' in merge and 'source_file_grid' in merge
    for relative in re.findall(r"fullfile\(repo_root,'results','raw','([^']+)','([^']+)'\)", merge):
        assert (root / 'results/raw' / relative[0] / relative[1]).is_file(), relative
    plot = (scripts / 'plot_entropy_maps_matlab.m').read_text(encoding='utf-8')
    assert plot.count("'processed', 'main_map'") == 3 and "'results', 'raw'" not in plot
    subprocess.run(['git','diff','--check'], cwd=root, check=True)
    bash = shutil.which('bash')
    if bash:
        for path in [*scripts.glob('*.sbatch'), *scripts.glob('*.sh')]:
            subprocess.run([bash, '-n', str(path)], check=True)
        print('PASS: Bash syntax checks')
    else:
        print('SKIP: Bash syntax checks (bash unavailable)')
    print('PASS: UTF-8, Python syntax, README commands, all main/convergence parameters, SLURM targets, test wiring, reconstruction inputs and output paths, git diff --check')


if __name__ == '__main__':
    main()
