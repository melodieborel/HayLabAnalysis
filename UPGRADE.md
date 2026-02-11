## Upgrade instructions (short)

1. Pull latest `main` after release:

```bash
git checkout main
git pull origin main
```

2. Run the smoke test to verify your environment:

```bash
pip install -e .[test]
pytest -q tests/test_smoke.py
```

3. Report any issues to the project channel and include the output of the smoke test.
