| Date | Step | Problem | Root Cause | Fix | Rerendered |
| --- | --- | --- | --- | --- | --- |
| 09/15/15 | Prealignment | 1,8 had slightly poor prealignment with 1,7 | 2/25 blockmatches failed to find the best match. In both cases, a dark patch in both images drove the correspondence, instead of the neurites. | Manually removed the bad matches | Yes |
| 09/16/15 | Montage | 1,40 failed to montage well | 1,40 is slightly out of focus, so cross correlation values were lower than the threshold in montage blockmatching | Reran with a lower threshold | Yes |
| 09/16/15 | Alignment | 1,2 failed locally in the interior aligning to 1,1 | There is a narrow scratch that runs through tiles 2 & 3  | Reran with a finer mesh | Yes |
| 09/17/15 | Prealignment | 1,33 failed to prealign well with 1,32 | 1,33 is translated by nearly 3000 px from 1,32 | Reran with a larger search radius | Yes |
| 09/19/15 | Montage | 1,103 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,120 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,130 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,151 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,153 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,154 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,155 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,165 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,166 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,167 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,52 montaged with a blank tile | Blank tile from the microscope was not removed in premontage & noisy image | Removed blank tile from premontage offset file | Yes |
| 09/19/15 | Montage | 1,54 montaged with a blank tile | Blank tile from the microscope was not removed in premontage & noisy image | Removed blank tile from premontage offset file | Yes |
| 09/20/15 | Prealignment | 1,54 failed to prealign due to lack of points | Large offset | Manually entered offset | Yes |
| 09/20/15 | Prealignment | 1,136 failed to prealign due to lack of points | Large offset | Manually entered offset | Yes |
| 09/20/15 | Prealignment | 1,139 failed to prealign due to lack of points | Large offset | Manually entered offset | Yes |
| 09/20/15 | Prealignment | 1,140 failed to prealign due to lack of points | Large offset | Manually entered offset | Yes |
| 09/22/15 | Montage | 1,56 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/28/15 | Prealignment | 1,119-121 has bad prealignment due to lack of points | Wide mesh | Reran with fine mesh | No |
| 09/28/15 | Prealignment | 1,129-132 has bad prealignment due to lack of points | Wide mesh | Reran with fine mesh | No |
| 09/28/15 | Prealignment | 1,141-143 has bad prealignment due to lack of points | Wide mesh + large offset | Reran with fine mesh / manually entered offset | No |
| 09/28/15 | Prealignment | 1,155-158 has bad prealignment due to lack of points | Large offset | Manually entered offset | No |
| 09/28/15 | Prealignment | 1,167-168 has bad prealignment due to lack of points | Large offset | Manually entered offset | No |
| 09/28/15 | Premontage | 2,3 has bad premontage due to an extra tile | NA | NA | No |
| 09/28/15 | Premontage | 2,52 has bad premontage due to an extra tile | NA | NA | No |
| 09/29/15 | Premontage | 3,22 has bad premontage due to an extra tile | NA | NA | No |
| 09/29/15 | Premontage | 3,90 has bad premontage | NA | NA | No |
| 09/29/15 | Premontage | 3,133 has bad premontage | NA | NA | No |
| 10/05/15 | Montage | 6,1 has bad montage | Missing tile | NA | No |
| 10/05/15 | Montage | 6,5 has bad montage | Missing tile | NA | No |
| 10/05/15 | Montage | 6,26 has bad montage | No overlap | NA | No |
| 10/05/15 | Montage | 6,53 has bad montage | No overlap | NA | No |
| 10/05/15 | Montage | 6,170 has bad montage | Little overlap | NA | No |
| 10/05/15 | Montage | 7,2 has bad montage | Missing tile | NA | No |
| 10/05/15 | Montage | 7,12 has bad montage | Missing tile | NA | No |
| 10/05/15 | Montage | 7,45 has bad montage | Little overlap (premontage) | NA | No |
| 10/05/15 | Montage | 7,50 has bad montage | Weird tile offsets | NA | No |
| 10/05/15 | Montage | 7,51 has bad montage | Missing tile | No | No |
| 10/05/15 | Montage | 7,63 has bad montage | Little overlap | No | No |
| 10/05/15 | Montage | 7,64 has bad montage | Little overlap | No | No |
| 10/05/15 | Montage | 7,70 has bad montage | Premontage failure | SKIP | No |
| 10/05/15 | Montage | 7,72 has bad montage | Premontage failure? | SKIP | No |
| 10/05/15 | Montage | 7,76 has bad montage | Little overlap | NA | No |
| 10/05/15 | Montage | 7,79 has bad montage | No overlap | NA | No |
| 10/05/15 | Montage | 7,87 has bad montage | Premontage far | NA | No |
| 10/05/15 | Montage | 7,99 has bad montage | Premontage far | NA | No |
| 10/05/15 | Montage | 7,102 has bad montage | Premontage far | NA | No |
| 10/05/15 | Montage | 7,105 has bad montage | Premontage far | NA | No |
| 10/05/15 | Montage | 7,117 has bad montage | No tiles | NA | No |
| 10/05/15 | Montage | 8,14 has bad montage | OOF | NA | No |
| 10/05/15 | Montage | 8,27 has bad montage | Small overlap | NA | No |
| 10/05/15 | Montage | 8,31-33 has bad montage | OOF | NA | No |
| 10/05/15 | Montage | 8,40 has bad montage | OOF | NA | No |
| 10/05/15 | Montage | 8,45 has bad montage | OOF | NA | No |
| 10/05/15 | Montage | 8,47 has bad montage | OOF | NA | No |
| 10/05/15 | Montage | 8,49 has bad montage | Small overlap | NA | No |
| 10/05/15 | Montage | 8,50 has bad montage | Small overlap | NA | No |
| 10/05/15 | Montage | 8,55 has bad montage | Small overlap | NA | No |
| 10/05/15 | Montage | 8,66 has bad montage | OOF | NA | No |
| 10/05/15 | Montage | 8,80 has bad montage | Missing tile | NA | No |
| 10/05/15 | Montage | 8,82 has bad montage | Missing tile | NA | No |
| 10/05/15 | Montage | 8,97 has bad montage | OOF | NA | No |
| 10/05/15 | Montage | 8,121 has bad montage | OOF | NA | No |
| 10/05/15 | Montage | 8,125 has bad montage | OOF | NA | No |
| 10/05/15 | Montage | 8,145 has bad montage | Missing tile | NA | No |
| 10/05/15 | Montage | 8,155 has bad montage | Missing tile | NA | No |
| 10/05/15 | Montage | 8,166 has bad montage | OOF | NA | No |





































