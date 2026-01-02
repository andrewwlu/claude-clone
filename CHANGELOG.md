# Changelog

All notable changes to this project will be documented in this file.

## [0.4.1] - 2025-12-28
### Fixed
- Overlap Tm calculation now uses nearest-neighbor method (was using simple GC method)
- Fixed edge case where fragments with identical 3' ends caused incorrect assembly order

### Changed
- Assembly Council now runs safety-expert with Sonnet instead of Haiku

## [0.4.0] - 2025-11-15
### Added
- Assembly Council: 10 parallel expert review system
- `/council-review` command for comprehensive design analysis
- Regulatory expert for compliance documentation

### Changed
- Migrated from local CSV to Google Sheets MCP integration

## [0.3.2] - 2025-09-03
### Fixed
- GenBank feature positions incorrect after rotation (issue #12)
- Homopolymer detection missed runs at fragment boundaries

## [0.3.0] - 2025-07-20
### Added
- MCP server integration for Google Sheets
- SessionStart hook for cache freshness check
- Audit logging for all file operations

### Deprecated
- `sync_from_csv()` - use MCP-based sync instead

## [0.2.1] - 2025-06-01
### Fixed
- Backbone cutting failed for enzymes with asymmetric overhangs

## [0.2.0] - 2025-05-10
### Added
- Automatic Twist order generation
- Fragment size validation against Twist constraints
- PostToolUse hook for GenBank validation

## [0.1.0] - 2025-03-15
### Added
- Initial release
- Basic Gibson assembly simulation
- GenBank file generation
- Fragment/Plasmid database sync from Google Sheets
