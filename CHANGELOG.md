# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.2] - 2019-09-12
## Fixed
- Change how lag indices used to hopefully prevent zero-division errors.

## [0.3.1] - 2017-07-18
## Changed
- Return np.array from anilagindices, instead of list of np.array.
- Refactor kriging.krige for readability.

## [0.3.0] - 2017-07-18
## Changed
- Replace zip(list(...)) with NumPy operations.

## [0.2.0] - 2019-07-18
### Fixed
- Update some import statements for Python3.
