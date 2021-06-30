# Elements

Elements currently exported by the Five

```@example
const LinearSolidElementQuad = LinearSolidElement{2,1,RefCube,Float64}
const LinearSolidElementTria = LinearSolidElement{2,1,RefTetrahedron,Float64}
const LinearSolidElementHexa = LinearSolidElement{3,1,RefCube,Float64}
const SolidElementQuad = SolidElement{2,1,RefCube,Float64}
const SolidElementHexa = SolidElement{3,1,RefCube,Float64}
const SolidElementTria = SolidElement{2,1,RefTetrahedron,Float64}
const SolidElementTetra = SolidElement{3,1,RefTetrahedron,Float64}
const RigidElement2d = RigidElement{2}
const RigidElement3d = RigidElement{3}
```

### Important functions

```@docs
Five.get_fields
Five.integrate_forcevector!
Five.integrate_forcevector_and_stiffnessmatrix!
Five.integrate_massmatrix!
```