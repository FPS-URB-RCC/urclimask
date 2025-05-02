# urclimask (Urban/Rural CLImate Mask defintion)

This repository provides a set of tools developed within the framework of the **CORDEX Flagship Pilot Study (FPS) on URBan environments and Regional Climate Change (URB-RCC)** to identify **urban** and surrounding **rural** areas from **climate model** outputs, tailored to specific cities of interest. The code includes functionalities to assess **Urban Heat Island (UHI)** intensity, with configurable parameters adaptable to different spatial resolutions and urban environments.

## Contents

| Directory | Contents |
| :-------- | :------- |
|  [urclimask]() | Python code to delimitate urban/rural mask and assees Urban Heat Island analysys.
|  [notebooks]() | Jupyter notebooks with examples on how to use the library for fiffenre RCMs and cities.
| [doc]() | Description of the model.

## Requirements

Scripts and (jupyter) notebooks are provided in [Python](https://www.python.org/) to ensure reproducibility and reusability of the results. The simplest way to match all these requirements is by using a dedicated [conda](https://docs.conda.io) environment, which can be easily installed by issuing:

```sh
conda create -n urclimask pip jupyter
conda activate urclimask
pip install urclimask
```

## Examples of use

Examples of use of the `urclimask` library are available in the form of [jupyter notebooks](). To run the examples follow the following steps:

1. Download the folder [notebooks]() from the github repository, or navigate to the folder should you have cloned the repo.
2. Open jupyter notebook of Jupyter Lab (type `jupyter notebook` or `jupyter lab`  in the terminal)
3. Open one of the tests available in the [notebooks]() folder with jupyter notebook  (e.g. [paris_across_CORDEX_resolutions.ipynb]())

## Errata and problem reporting

To report an issue with the library, please fill a GitHub issue.

## License
Copyright 2023, European Union.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
