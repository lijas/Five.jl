export LinearSolidElementHexa, LinearSolidElementQuad, LinearSolidElementTria

#Definition of all elements
const LinearSolidElementQuad = LinearSolidElement{2,1,RefCube,Float64}
const LinearSolidElementTria = LinearSolidElement{2,1,RefTetrahedron,Float64}
const LinearSolidElementHexa = LinearSolidElement{3,1,RefCube,Float64}
const SolidElementQuad = SolidElement{2,1,RefCube,Float64}
const SolidElementHexa = SolidElement{3,1,RefCube,Float64}
const SolidElementTria = SolidElement{2,1,RefTetrahedron,Float64}
const SolidElementTetra = SolidElement{3,1,RefTetrahedron,Float64}
const RigidElement2d = RigidElement{2}
const RigidElement3d = RigidElement{3}
#const SpringElement2d = SpringElement{2}
#const MyBeamElement = BeamElement{Float64}
const ShellElementQuad = ShellElement{Float64}
const BLTShell = BelytschkoLinTsayShellElement{Float64}

#const ALL_ELEMENTS = [SpringElement2d, LinearSolidElementQuad, LinearSolidElementTria, SolidElementQuad, SolidElementTria, RigidElement2d]