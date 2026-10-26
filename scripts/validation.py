"""Shared execution/provenance utilities; no generated performance observations."""
import hashlib,json,os,shutil,subprocess,sys,threading,time
from pathlib import Path

def digest(path):
    with open(path, 'rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def toolchain_hashes(directory):
    """Bind the launcher, actual reader and nested runtime resources together."""
    return {str(path.relative_to(directory)):digest(path) for path in sorted(directory.rglob('*'))
            if path.is_file() and '__pycache__' not in path.parts and path.suffix!='.pyc'}

class Runner:
    def __init__(self, root, binpath, *, resume=False):
        self.root, self.rows = root, []
        # A run must not silently pick up executables rebuilt during development.
        self.binpath = root / 'toolchain'
        if resume:
            if binpath.resolve() != self.binpath.resolve() or not self.binpath.is_dir():
                raise RuntimeError('resume requires the existing frozen toolchain')
            self.rows = json.loads((root / 'commands.json').read_text())
        else:
            self.binpath.mkdir()
            for name in ('FastGA-baseline','FastGA-shared','FastGA-shared-core','GIXmake','GIXmake-direct','GIXmake-marks','pfp-emit','joint-copy','FAtoGDB','ALNtoPAF','ALNtoPSL','gixshare','build.json'):
                source = binpath / name
                if source.is_file():
                    shutil.copy2(source,self.binpath/name)
            for name in ('fastga-runtime','pfp-tools'):
                if (binpath/name).is_dir():
                    shutil.copytree(binpath/name,self.binpath/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        self.env = dict(os.environ, PATH=str(self.binpath) + os.pathsep + os.environ['PATH'])

    def run(self, name, args, env=None, expected_error=None, measure_tree=False, temporary_dir=None):
        if measure_tree:
            import psutil
        stdout, stderr = self.root / (name + '.stdout'), self.root / (name + '.stderr')
        current = dict(self.env)
        current.update(env or {})
        start = time.monotonic()
        with stdout.open('wb') as out, stderr.open('wb') as err:
            proc = subprocess.Popen(list(map(str, args)), stdout=out, stderr=err,
                                    cwd=self.root, env=current)
            peak = [0];scratch_peak=[0,0];scratch_samples=[0];stop = threading.Event()
            def sample():
                process=None
                if measure_tree:
                    try: process = psutil.Process(proc.pid)
                    except psutil.NoSuchProcess: pass
                while not stop.is_set():
                    if process is not None:
                        try:
                            processes = [process] + process.children(recursive=True)
                            total = 0
                            for child in processes:
                                try: total += child.memory_info().rss
                                except psutil.NoSuchProcess: pass
                            peak[0] = max(peak[0],total)
                        except psutil.NoSuchProcess: pass
                    if temporary_dir is not None:
                        logical=allocated=0;seen=set()
                        for path in Path(temporary_dir).rglob('*'):
                            try:
                                if path.is_symlink() or not path.is_file():continue
                                st=path.stat();key=(st.st_dev,st.st_ino)
                                if key in seen:continue
                                seen.add(key);logical+=st.st_size;allocated+=getattr(st,'st_blocks',0)*512
                            except FileNotFoundError:continue
                        scratch_peak[0]=max(scratch_peak[0],logical);scratch_peak[1]=max(scratch_peak[1],allocated);scratch_samples[0]+=1
                    stop.wait(.5)
            monitor = threading.Thread(target=sample,daemon=True) if measure_tree or temporary_dir is not None else None
            if monitor: monitor.start()
            try:
                _, status, usage = os.wait4(proc.pid, 0)
                proc.returncode = os.waitstatus_to_exitcode(status)
            finally:
                stop.set()
                if monitor: monitor.join()
        row = dict(name=name, command=list(map(str, args)), exit_code=proc.returncode,
                   wall_seconds=time.monotonic() - start,
                   peak_rss_bytes=usage.ru_maxrss * (1 if sys.platform == 'darwin' else 1024),
                   stdout=str(stdout), stderr=str(stderr))
        if measure_tree:
            row.update(sampled_process_tree_peak_rss_bytes=peak[0] or None,rss_sample_interval_seconds=.5,
                       rss_note='Sum of resident sets; shared pages may be counted more than once.')
        if temporary_dir is not None:row.update(sampled_peak_scratch_logical_bytes=scratch_peak[0],sampled_peak_scratch_allocated_bytes=scratch_peak[1],scratch_sample_interval_seconds=.5,scratch_samples=scratch_samples[0],scratch_scope=str(temporary_dir))
        self.rows.append(row)
        (self.root / 'commands.json').write_text(json.dumps(self.rows, indent=2) + '\n')
        if expected_error:
            if proc.returncode == 0 or expected_error not in stderr.read_text():
                raise RuntimeError(f'{name}: expected labelled negative test rejection was absent')
            row['EXPECTED_TEST_REJECTION'] = expected_error
        elif proc.returncode:
            raise RuntimeError(f"Program {args[0]} failed with: {stderr.read_text()[-5000:]}")
        return stdout

def gix_files(root):
    return sorted([*root.glob('*.gix'), *root.glob('.*.ktab.*')])
