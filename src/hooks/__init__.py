"""
Hook scripts for Claude Code integration.

This module contains scripts that are executed by Claude Code hooks:
- session_init: Run on SessionStart
- session_archive: Run on SessionEnd
- validate_genbank: Run on PostToolUse(Write) for .gb files
- validate_council_output: Run on SubagentStop
- audit_log: Run on PostToolUse for logging
- pre_write_validator: Run on PreToolUse(Write)
"""
