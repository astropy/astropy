# NOTE: this hook should be added to
# https://github.com/pyinstaller/pyinstaller-hooks-contrib
# once that repository is ready for pull requests
from PyInstaller.utils.hooks import collect_data_files

datas = collect_data_files('skyfield')
