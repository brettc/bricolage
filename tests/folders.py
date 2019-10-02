import pathlib
import sys

# Add the parent here
here_dir = pathlib.Path(__file__).parent.absolute()
parent_dir = here_dir.parent
data_dir = here_dir / "data"
sys.path.append(str(parent_dir))
