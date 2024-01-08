Xray tools

add this repo to your system $PATH and $PYTHONPATH

convert volumes to zarr files

vol_to_zarr.py <fn.vol> z x y chunk

e.g.

vol_to_zarr.py volume_1.vol 200 3000 3000 100

vol_to_zarr.py volume_2.vol 200 3000 3000 100

run resolution estimation on the volumes

resolution_measure.py <zarr file1> <zarr file2> <ncores> <cube size>[ -ps <pixel_size=50> --snrt <snrt value=.2071>]

resolution_measure.py jaspersLegCryo_r1_50nm_rec_cone_01799_10001800.zarr jaspersLegCryo_r1_50nm_rec_cone_12000_10001800.zarr 8 200 --snrt 0.143
