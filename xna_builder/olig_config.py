from pathlib import Path
from typing import Any
from pydantic import BaseModel
import json


class OligoConfig(BaseModel):
    """
    Configuration model for building a single oligonucleotide structure.
    """

    name: str
    seq_base: str
    seq_backbone: list[str]
    save_dir: Path
    pnab_config: Path
    pnab_kwargs: dict[str, Any]

    @classmethod
    def from_json(cls, config_path: Path) -> list["OligoConfig"]:
        """
        Loads configurations from a JSON file and resolves all relative paths
        to be absolute paths relative to the config file's directory.

        Args:
            config_path (Path): The path to the input JSON configuration file.

        Returns:
            List[OligoConfig]: A list of validated OligoConfig objects.
        """
        config_dir = config_path.resolve().parent

        validated_configs = []

        with open(config_path, "r") as f:
            raw_configs = json.load(f)

        if not isinstance(raw_configs, list):
            raw_configs = [raw_configs]

        for i, raw_config in enumerate(raw_configs):

            if "save_dir" in raw_config:
                raw_config["save_dir"] = config_dir / raw_config["save_dir"]

            if "pnab_config" in raw_config:
                raw_config["pnab_config"] = config_dir / raw_config["pnab_config"]

            config = cls.model_validate(raw_config)
            validated_configs.append(config)

        return validated_configs
