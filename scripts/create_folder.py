import argparse
from pathlib import Path


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Checks whether a directory exists and creates it if necessary."
    )

    parser.add_argument(
        "directory",
        type=Path,
        help="Path to the directory to be created.",
    )

    return parser.parse_args()


def main():
    args = parse_arguments()

    directory = args.directory

    if not directory.is_dir():
        directory.mkdir(parents=True, exist_ok=True)
        print(f"Directory created: {directory}")
    else:
        print(f"Directory already exists: {directory}")


if __name__ == "__main__":
    main()
