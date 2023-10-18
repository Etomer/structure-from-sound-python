# Structure from Sound – python

## Setup
### Requirements
Install requirements with 

```
pip install -r requirements.txt
```

### Download Data
- Download data from [here](https://vision.maths.lth.se/erik_test/) and put the files in the `data` folder


### Using GUI-plot

In order to use `gui_plot.py` first run `genereate_detection_data.py` (takes ~10 min). Then run `gui_plot.py`

### Structure sketch
![structure sketch](./images/structure_sketch.png)

#### Data Formats

- Audio is stored as `.wav` files where each experiment is stored in a serperate folder with microphones i recording stored as `Track i.wav`
- Detections are stored as `<experiment name>.npy` with the stored array being a 2D matrix where each row is a detection and each column has the following information `<mic 1 index>, <mic 2 index>, <time t from start of recording>, <tdoa-detection (i.e. d_ijt = |s(t) - r_2| - |s(t) - r_1|)>`

