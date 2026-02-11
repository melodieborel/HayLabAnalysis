## Release notes for migration/refactor-core (draft)

Summary
- Move core utilities into `HayLabAnalysis.core` and keep adapter shims for compatibility.

User impact
- No changes required for normal use. Public imports like `from HayLabAnalysis import tools` continue to work.
- A deprecation warning will be emitted on import; this will be removed in a future minor release.

Rollback
- Revert the merge commit if an urgent rollback is needed.
