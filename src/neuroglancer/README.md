# Alembic Inspect
Inspect and edit correspondences from Match objects.

## Requirements

### Python
* python2.7

### Web app
* [Seung Lab fork of Neuroglancer](https://neuromancer-seung-import.appspot.com/).

## Getting started
1. `pip install [PACKAGEDIR]/requirements.txt`  

1. Within python:
```
using StateServer
pc = PointsController(9999) # your favorite port
```  
You're now broadcasting a state server at this address: `https://localhost:9999`.  
1. This state server is not secure, so your browser requires permission to access it. Open your browser and go to your state server's address, `https://localhost:9999`, and click through prompts to grant permission.
1. Still in your browser, load up your dataset of interest in neuroglancer, including the URL to the state server you've created in its state. For example:
`https://neuromancer-seung-import.appspot.com/#!{'layers':{'image':{'type':'image'_'source':'precomputed://gs://neuroglancer/pinky40_v11/image'}}_'navigation':{'pose':{'position':{'voxelSize':[4_4_40]_'voxelCoordinates':[60126.37109375_10611.818359375_782]}}_'zoomFactor':44.34634372224735}_'layout':'xy'_'stateServer':'https://localhost:9999'}`  
See how your state server address is passed to the 'stateServer' key at the end of the URL.
1. Back to Julia, you now have access to the points layer created in neuroglancer.
```
pc[:get]() # returns []
pc[:set]([])
```


## Credits
Thanks to [the original Neuroglancer project team](https://github.com/google/neuroglancer)!