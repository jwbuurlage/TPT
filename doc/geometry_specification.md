# Acquisition geometry specification

We specify acquisition geometries using a [TOML](https://github.com/toml-lang/toml/blob/master/versions/en/toml-v0.4.0.md) file. An example is given here:

```toml
specifies = "geometry"
kind = "parallel"
dimension = 3

[volume]
origin = [0.0, 0.0, 0.0]
lengths = [10.0, 10.0, 10.0]
voxels = [100, 100, 100]

[parameters]
source-distance = 5.0
detector-distance = 5.0
projection-count = 10
# angles: [0, 0.5, 1.0, ...]
# ...
```

The geometry specs have the following structure:

- First, we specify that the file is a geometry specification:
    ```toml
    specifies = geometry
    ```
- Then follows the *class* of the geometry as:
    ```toml
    class = name
    ```
    where name is one of the following:
    * `parallel`
    * `cone-beam`
    * `laminography`
    * `tomosynthesis`
    * `...`
- Next the problem dimension is given as:
    ```toml
    dimension = 3
    ```
- Next a *volume* table is given, this table has *three* entries:
    * `origin` should be an array of `D` floating point numbers. This is the left-most point of the volume (in each axis).
    * `lengths` should be an array of `D` floating point numbers, representing the physical dimensions.
    * `voxels` should be an array of `D` integers, the number of voxels along each axis.
- Next a table of parameters are given. These depend on the geometry and are outlined below

## Geometry parameters

### Parallel

- `source-distance`
- `detector-distance`
- (optional) `projection-count`
- (optional) `angles`
