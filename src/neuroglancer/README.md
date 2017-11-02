# Inspect
Inspect and edit correspondences from Match objects.

## Requirements

### Python
* python2.7

### Web app
* [Seung Lab fork of Neuroglancer](https://neuromancer-seung-import.appspot.com/).

## Getting started
1. `pip install -r [PACKAGEDIR]/requirements.txt`  

1. Open python
```python
```
1. From within the REPL
```from match_inspect import *
ms = MeshSetInspect([PATH TO MESH CSV FILES], 9999) # your favorite port
``` 
You should see a statement 'IOLoop starting'. You're now broadcasting a state server at this address: `https://localhost:9999`.  
1. This state server is not secure, so your browser requires permission to access it. Open your browser and go to your state server's address, `https://localhost:9999`, and click through prompts to grant permission.
1. Still in your browser, load up your dataset of interest in neuroglancer, including the URL to the state server you've created in its state. For example:
`https://neuromancer-seung-import.appspot.com/#!{'layers':{'image':{'type':'image'_'source':'precomputed://gs://neuroglancer/pinky40_v11/image'}}_'navigation':{'pose':{'position':{'voxelSize':[4_4_40]_'voxelCoordinates':[60126.37109375_10611.818359375_782]}}_'zoomFactor':44.34634372224735}_'layout':'xy'_'stateServer':'https://localhost:9999'}`  
See how your state server address is passed to the 'stateServer' key at the end of the URL.
1. Make a change to the state of your neuroglancer app, e.g. pan in the image.
1. Back in python, you should see a statement 'state initialized', along with a print out of the JSON describing the state.
1. You can now interact with your matches as follows.
```ms.set_id(0);		# set to inspec the first match csv in the directory
m = ms.inspect()	# create a MatchInspect object that will display point-pairs
```
1. Within neuroglancer, you should now see point-pairs. You can `CTRL+left-click` to remove them. (You can also use `CTRL+left-click to add them, but the MatchInspect object will ignore new points.)
1. In python, you can flip the labels so that all rejected points become accepted, and vice versa.
```m.flip()
```
1. When the matches you have in the neuroglancer app are correct, you can save them back to their csv with:
```ms.save()
```
1. You can move on to the next match file with
```m = ms.get_next()
```

## Credits
Thanks to [the original Neuroglancer project team](https://github.com/google/neuroglancer)!