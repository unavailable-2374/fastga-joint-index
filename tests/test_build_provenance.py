#!/usr/bin/env python3
"""Check actual native build manifests and isolated checksum/tamper rejections."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from validation import digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bin',type=Path,required=True)
    p.add_argument('--report',type=Path,required=True)
    a=p.parse_args();build=a.bin.resolve()
    evidence=dict(status='RUNNING',kind='BUILD VERIFICATION; no biological data',checks=[])
    def checked(label):evidence['checks'].append(dict(check=label,status='PASS'))
    try:
        manifest=json.loads((build/'build.json').read_text())
        assert manifest['build_system']=='CMake' and manifest['build_python_required'] is False
        assert manifest['provenance_generator_sha256']==digest(ROOT/'tools/build_manifest.cpp')
        for compiler,version in manifest['compiler_versions'].items():
            result=subprocess.run([compiler,'--version'],capture_output=True,text=True)
            if result.returncode:raise RuntimeError('Program '+compiler+' failed with: '+result.stderr)
            assert result.stdout==version
        checked('actual compiler versions and native provenance writer')
        for name,sha in manifest['frontend_sha256'].items():assert digest(build/name)==sha
        for name,sha in manifest['adapter_sha256'].items():assert digest(ROOT/'src'/name)==sha
        for name,sha in manifest['source_sha256'].items():assert digest(build/'upstream/fastga'/name)==sha
        pfp=json.loads((build/'pfp-tools/build.json').read_text())
        for name,sha in pfp['executables'].items():assert digest(build/'pfp-tools'/name)==sha
        checked('all reported source, frontend and PFP hashes match actual bytes')
        with tempfile.TemporaryDirectory(prefix='fastga-vendor-QC-') as directory:
            root=Path(directory)
            shutil.copytree(ROOT/'vendor',root/'vendor')
            shutil.copytree(ROOT/'patches',root/'patches')
            shutil.copytree(ROOT/'cmake',root/'cmake')
            (root/'CMakeLists.txt').write_text('cmake_minimum_required(VERSION 3.20)\nproject(VendorQC NONE)\ninclude(cmake/PrepareVendor.cmake)\n')
            command=['cmake','-S',str(root),'-B',str(root/'build')]
            def configure(error=None):
                result=subprocess.run(command,capture_output=True,text=True)
                if error:
                    assert result.returncode and error in result.stderr
                    checked('EXPECTED TEST REJECTION: '+error)
                elif result.returncode:
                    raise RuntimeError('Program CMake failed with: '+result.stderr+result.stdout)
                return result
            configure()
            files=[p for p in (root/'build/upstream').rglob('*') if p.is_file()]
            mtimes={p:p.stat().st_mtime_ns for p in files}
            result=configure()
            assert 'reused=TRUE' in result.stdout
            assert mtimes=={p:p.stat().st_mtime_ns for p in files}
            checked('reconfigure verifies and reuses unchanged patched sources')
            prepared=root/'build/upstream/fastga/FastGA.c';expected=digest(prepared)
            prepared.write_bytes(prepared.read_bytes()+b'\n/* TEST tampering, isolated copy */\n')
            result=subprocess.run(['cmake','--build',str(root/'build')],capture_output=True,text=True)
            if result.returncode:raise RuntimeError('Program CMake failed with: '+result.stderr+result.stdout)
            assert digest(prepared)==expected
            checked('build automatically rechecks and restores tampered build-local source')
            source=root/'vendor/fastga/FastGA.c';original=source.read_bytes()
            source.write_bytes(original+b'\n/* TEST tampering, isolated copy */\n')
            configure('vendored source checksum mismatch: fastga/FastGA.c')
            source.write_bytes(original);configure()
            checked('valid pinned source recovery after expected rejection')
        evidence['status']='PASS'
    except BaseException as error:
        evidence.update(status='FAILED',error=str(error));raise
    finally:a.report.write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps(dict(status='PASS',checks=len(evidence['checks']))))


if __name__=='__main__':main()
