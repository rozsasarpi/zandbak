# zandbak

[![Documentation Status](https://readthedocs.org/projects/zandbak/badge/?version=latest)](https://zandbak.readthedocs.io/en/latest/?badge=latest)
[![Continuous integration](https://github.com/rozsasarpi/zandbak/actions/workflows/push.yaml/badge.svg)](https://github.com/rozsasarpi/zandbak/actions)
[![PyPI version](https://img.shields.io/pypi/v/zandbak)](https://pypi.org/project/zandbak/)
[![coverage](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/rozsasarpi/da9e3419b54a0daf6fe07b934f37f837/raw/zandbak_main_coverage.json)](https://en.wikipedia.org/wiki/Code_coverage)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/rozsasarpi/zandbak.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/rozsasarpi/zandbak/context:python)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)


A sandbox python repository for testing GitHub features, CI pipelines, and python
packaging.

Some structural mechanics examples are used as dummy code:
 * computing the responses of linear elastic beams.


### Documentation

[Documentation](https://zandbak.readthedocs.io/en/latest/?badge=latest) based on
the latest `main`.


### Code coverage

* Solely based on using GitHub Actions, based on [this guide](https://dev.to/thejaredwilcurt/coverage-badge-with-github-actions-finally-59fa).
* Upon each push to any of the branches a json file with coverage information is
  added to [a dedicated gist](https://gist.github.com/rozsasarpi/da9e3419b54a0daf6fe07b934f37f837).
* The content of the json file is used to generate a badge using https://shields.io/.
* This badge is offered as part of the pull request template.
* If you use the same branch name multiple times then the old coverage badges will be
  overwritten, and in turn the old pull request badges will be incorrect (it is deemed
  to be still acceptable since the coverage of closed pull requests are not expected to be
  checked, but if that is needed then the correct coverage information can be obtained
  from the pipeline prints).
