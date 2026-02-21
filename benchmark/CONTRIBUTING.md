# Contributing to MagNav TL Benchmark Framework

We welcome contributions from the community! This document provides guidelines for contributing to the MagNav TL Benchmark Framework.

## Table of Contents
- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Workflow](#development-workflow)
- [Submitting Changes](#submitting-changes)
- [Reporting Issues](#reporting-issues)
- [Documentation](#documentation)

## Code of Conduct

This project follows a code of conduct to ensure a welcoming environment for all contributors. By participating, you agree to:
- Be respectful and inclusive
- Focus on constructive feedback
- Accept responsibility for mistakes
- Show empathy towards other contributors

## Getting Started

### Prerequisites
- Julia 1.10 or later
- Git
- Basic understanding of magnetic navigation and Tolles-Lawson compensation

### Development Setup

1. **Fork the repository** on GitHub
2. **Clone your fork**:
   ```bash
   git clone https://github.com/yourusername/magnav-benchmark.git
   cd magnav-benchmark
   ```

3. **Set up the development environment**:
   ```bash
   julia --project=. -e "using Pkg; Pkg.instantiate()"
   ```

4. **Run tests** to ensure everything works:
   ```bash
   julia --project=. runtests.jl
   ```

5. **Create a feature branch**:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## Development Workflow

### Code Standards

- **Julia Style**: Follow the [Julia Style Guide](https://docs.julialang.org/en/v1/manual/style-guide/)
- **Documentation**: Document all public functions with docstrings
- **Type Hints**: Use appropriate type annotations
- **Error Handling**: Provide meaningful error messages
- **Performance**: Consider computational efficiency for large datasets

### Testing

- Add unit tests for new functionality in the `test/` directory
- Ensure all tests pass before submitting changes
- Test on multiple Julia versions when possible
- Include integration tests for complex workflows

### File Organization

```
magnav-benchmark/
├── Core algorithms in root directory
│   ├── custom_tl.jl          # TL model implementations
│   ├── data_loader.jl        # Data loading utilities
│   └── nav_metrics.jl        # Performance metrics
├── Test files
│   ├── test_*.jl             # Unit tests
│   └── runtests.jl           # Test runner
└── Documentation
    ├── README.md             # Main documentation
    └── docs/                 # Additional docs
```

## Submitting Changes

### Pull Request Process

1. **Update documentation** for any changed functionality
2. **Add tests** for new features
3. **Ensure CI passes** (if available)
4. **Update CHANGELOG.md** with your changes
5. **Submit a pull request** with a clear description

### Commit Messages

Use clear, descriptive commit messages:
```
feat: add new TL model variant with cross-validation
fix: correct heading correlation calculation in nav_metrics.jl
docs: update installation instructions for Julia 1.11
```

### Pull Request Template

Please include:
- **Description**: What changes were made and why
- **Testing**: How the changes were tested
- **Breaking Changes**: Any breaking changes and migration guide
- **Screenshots**: For UI/visualization changes

## Reporting Issues

### Bug Reports

When reporting bugs, please include:
- **Julia version** and platform
- **Minimal reproducible example**
- **Expected vs actual behavior**
- **Error messages and stack traces**
- **Data sample** (if applicable, anonymized)

### Feature Requests

For feature requests, please:
- **Describe the problem** you're trying to solve
- **Explain your proposed solution**
- **Consider alternative approaches**
- **Discuss potential impact** on existing functionality

## Documentation

### Code Documentation

All public functions should have docstrings:
```julia
"""
    fit_custom_tl(xyz, cal_indices; λ=0.0)

Fit custom 9-term Tolles-Lawson model to calibration data.

# Arguments
- `xyz`: XYZ20 struct with flight data
- `cal_indices`: Indices for calibration segment
- `λ`: Ridge regularization parameter (default: 0.0)

# Returns
- Coefficient vector for the 9-term model

# Examples
```julia
coef = fit_custom_tl(xyz, 1:6000)
```
"""
function fit_custom_tl(xyz, cal_indices; λ=0.0)
    # implementation...
end
```

### README Updates

When adding new features:
- Update the main README.md with usage examples
- Add to the API reference section
- Update installation instructions if needed
- Include performance benchmarks for new algorithms

## Recognition

Contributors will be recognized:
- In CHANGELOG.md for all contributions
- In the README acknowledgments section for significant contributions
- As co-authors on related publications

Thank you for contributing to the MagNav TL Benchmark Framework! 🚀