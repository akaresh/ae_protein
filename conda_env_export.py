
"""
Export a Conda environment with --from-history, but also append
Pip-installed dependencies
Exports only manually-installed dependencies, excluding build versions, but
including Pip-installed dependencies.
Lots of issues requesting this functionality in the Conda issue tracker, no
sign of progress (as of March 2020).
TODO (?): support command-line flags -n and -p
"""
import re
import subprocess
import sys

import yaml


def export_env(history_only=False, include_builds=False):
    """ Capture `conda env export` output """
    cmd = ['conda', 'env', 'export']
    if history_only:
        cmd.append('--from-history')
        if include_builds:
            raise ValueError(
            	'Cannot include build versions with "from history" mode')
    if not include_builds:
        cmd.append('--no-builds')
    cp = subprocess.run(cmd, stdout=subprocess.PIPE)
    try:
        cp.check_returncode()
    except:
        raise
    else:
        return yaml.safe_load(cp.stdout)


def _is_history_dep(d, history_deps):
    if not isinstance(d, str):
        return False
    d_prefix = re.sub(r'=.*', '', d)
    #print(d_prefix)
    if d_prefix in history_deps: return True
    else: return False
    #return d_prefix in history_deps


def _get_pip_deps(full_deps):
    for dep in full_deps:
        if isinstance(dep, dict) and 'pip' in dep:
            return dep


def _combine_env_data(env_data_full, env_data_hist):
    deps_full = env_data_full['dependencies']
    deps_hist = env_data_hist['dependencies']
    deps = [dep for dep in deps_full if _is_history_dep(dep, deps_hist)]
    #print(deps)
    #sys.exit()
    pip_deps = _get_pip_deps(deps_full)
    
    env_data = {}
    env_data['channels'] = env_data_full['channels']
    env_data['dependencies'] = deps
    env_data['dependencies'].append(pip_deps)

    return env_data


def main(name=None, python=None):
    env_data_full = export_env()
    env_data_hist = export_env(history_only=True)
    env_data = _combine_env_data(env_data_full, env_data_hist)
    if name is not None: env_data['name'] = name
    if python is not None:
    	assert('=' in python)
    	env_data['dependencies'].append(python)
    
    yaml.dump(env_data, sys.stdout)
#     print(
#     """
#     Warning: this output might contain packages installed from
#     non-public sources, e.g. a Git repository.
#     You should review and test the output to make sure it works with 
#     `conda env create -f`, 'and make changes as required.
#     'For example, `conda-env-export` itself is not currently uploaded to 
#     PyPI, and it must be removed from the output file, or else 
#     `conda create -f` will fail.', file=sys.stderr
#     """)


if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(
		description='Export conda env yaml file with pip dependencies')
	parser.add_argument(
		'--name', '-n', required=False, type=str,
		metavar='<str>', help='name for conda environment')
	parser.add_argument(
		'--python', '-p', required=False, type=str,
		metavar='<str>', default='python=3.8', help='python version')	
	
	arg = parser.parse_args()
	main(name=arg.name, python=arg.python)