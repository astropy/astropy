from sherpa import models






conda create -n astropy_with_sherpabridge2 python=2.7  numpy cython jinja2
git clone https://github.com/nocturnalastro/astropy.git astropy_with_sherpabridge
conda install -c https://conda.binstar.org/sherpa sherpa
source-bash source activate astropy_with_sherpabridge2
cd astropy_with_sherpabridge
git checkout sherpa_bridge
python setup.py install
python -c "from astropy.modeling.astro_sherpa import SherpaFitter"
python setup.py test