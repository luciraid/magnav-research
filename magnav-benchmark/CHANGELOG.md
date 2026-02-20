# Changelog

All notable changes to the MagNav TL Benchmark Framework will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2024-02-20

### Added
- **Complete benchmarking framework** for custom 9-term vs official MagNav TL compensation
- **Automatic data loading** via Julia artifacts (SGL 2020 flight data + Ottawa maps)
- **Enhanced plotting suite** with publication-quality visualizations and statistical annotations
- **Interactive Pluto notebook** support for exploratory analysis
- **Comprehensive documentation** with quick-start guide and API reference
- **Unit testing framework** with validation scripts
- **Multi-segment cross-validation** for robust performance assessment
- **Navigation performance evaluation** with Extended Kalman Filter integration
- **Residual analysis tools** for TL compensation quality assessment

### Features
- **Data Processing**: 108K+ flight samples with 10Hz sampling rate
- **TL Models**: Custom 9-term implementation vs MagNav official implementation
- **Performance Metrics**: DRMS, CEP, heading correlation, RMS improvement
- **Visualization**: Error time series, position scatter plots, heading correlation analysis
- **Reproducibility**: Complete parameter logging and artifact-based data management

### Technical Details
- **Julia Version**: 1.10+ compatibility
- **Dependencies**: MagNav.jl, Plots.jl, HDF5.jl, LinearAlgebra.jl
- **Data Sources**: SGL 2020 training dataset, Ottawa magnetic anomaly maps
- **Output Formats**: PNG/PDF plots, CSV metrics, timestamped result directories

### Performance Results
- **TL Compensation**: 40%+ variance reduction in calibration segments
- **Magnetic Field Improvement**: 25-35% RMS reduction in evaluation segments
- **Navigation Enhancement**: Quantifiable position accuracy improvements

### Documentation
- **README.md**: Complete installation and usage guide
- **API Reference**: All public functions documented
- **Troubleshooting Guide**: Common issues and solutions
- **Contributing Guidelines**: Development workflow and standards

## [0.1.0] - 2024-01-15

### Added
- Initial project structure and core components
- Basic TL model implementations
- Data loading infrastructure
- Unit test framework
- Preliminary documentation

### Changed
- Migrated from manual data management to Julia artifacts
- Updated API compatibility for MagNav.jl v0.2+
- Improved error handling and validation

### Fixed
- Memory issues with large flight datasets
- API compatibility problems with XYZ20 struct
- Unicode encoding issues in source files

---

## Types of changes
- `Added` for new features
- `Changed` for changes in existing functionality
- `Deprecated` for soon-to-be removed features
- `Removed` for now removed features
- `Fixed` for any bug fixes
- `Security` in case of vulnerabilities

## Versioning Policy
This project follows [Semantic Versioning](https://semver.org/):
- **MAJOR** version for incompatible API changes
- **MINOR** version for backwards-compatible functionality additions
- **PATCH** version for backwards-compatible bug fixes

---

*For older changes, see the git history.*