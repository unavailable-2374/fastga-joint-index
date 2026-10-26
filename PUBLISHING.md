# Repository and publishing

The canonical repository is
<https://github.com/unavailable-2374/fastga-joint-index>. The default development
branch is `main`.

```sh
git clone https://github.com/unavailable-2374/fastga-joint-index.git
cd fastga-joint-index
python3 -m pip install -r requirements-benchmark.txt
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSJ_PYTHON_TESTS=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
```

For subsequent changes, commit the intended source files and push the branch:

```sh
git push -u origin main
```

Before public distribution, choose a license for the new integration code.
Third-party license texts and notices are already retained. The repository does
not claim that an integration-wide license has been selected.

Generated build folders, genotype sequences, indexes, large alignment outputs
and private cluster details are excluded. `benchmarks/` contains explicitly
scoped real summaries and separate standalone test evidence. Human speed and
merged self-alignment measurements remain incomplete because those experiments
were stopped by request.

The Linux/macOS GitHub Actions workflow is included with pinned action commits.
Check the repository's [Actions page](https://github.com/unavailable-2374/fastga-joint-index/actions)
for its current results. Verification from the original local extraction is
recorded in `benchmarks/standalone-verification.json`; its CI field describes
the state at that local verification, before the initial upload.
