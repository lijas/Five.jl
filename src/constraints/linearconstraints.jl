
struct LinearConstraint{T}
	dofs::Vector{Int}
	coeffs::Vector{T}
	g::T
end

mutable struct LinearConstraints{dim,T}
	dh::DofHandler{dim,T}
	constraints::Vector{LinearConstraint{T}}

	ddofs::Vector{Int} #dependent dofs
	idofs::Vector{Int} #independnt dofs
	rdofs::Vector{Int} #the Rest of the dofs

	nconstraints::Ref{Int}

	C::Matrix{T}
	G::Vector{T}
end

function LinearConstraints{dim,T}(dh::DofHandler) where dim where T
	return LinearConstraints{dim,T}(dh, LinearConstraint[], Int[],Int[],Int[], 0, zeros(T,1,1), zeros(T,1))
end

function add_lc!(lc, coeffs, g, dofs)
	push!(lc.constraints, LinearConstraint(dofs,coeffs,g))
	lc.nconstraints[] += 1
end

function close_lc!(lc)

	I = Int[]
	J = Int[]
	V = Float64[]

	#Get all dofs that are used in the constraints
	_dofs = Int[]
	for c in lc.constraints
		for d in c.dofs
			push!(_dofs, d)
		end
	end


	constraint_dofs = unique(_dofs)
	dof_mapper = Dict{Int,Int}()

	#Create dof_mapper global_dof -> idx
	for (i,d) in enumerate(constraint_dofs)
		dof_mapper[d] = i
	end

	G = zeros(lc.nconstraints[])
	for (i,c) in enumerate(lc.constraints)		
		for j in 1:length(c.dofs)
			push!(I, i)
			push!(J, dof_mapper[c.dofs[j]])
			push!(V, c.coeffs[j])
		end
		G[i] = c.g
	end

	C = sparse(I, J, V)

	L,U,dp = lu(Matrix(C))

	#dp = [2,3]

	ip = setdiff(1:length(constraint_dofs),dp)

	lc.ddofs = constraint_dofs[dp]
	lc.idofs = setdiff(constraint_dofs, lc.ddofs)
	lc.rdofs = setdiff(1:ndofs(lc.dh), constraint_dofs)

	#Calculate lc.C, 
	#dependent_dofs = lc.C * independent_dofs
	lc.C = -inv(Matrix(C[:,dp]))*C[:,ip]
	lc.G = -inv(Matrix(C[:,dp]))*G[dp]

end

function apply_lc!(K, f, lc)
	C = lc.C
	G = lc.G
	idofs = lc.idofs; rdofs = lc.rdofs; ddofs = lc.ddofs

	f[rdofs] -= K[rdofs, ddofs]*G
	f[idofs] -= K[idofs, ddofs]*G + C'*f[ddofs] - C'*K[ddofs,ddofs]*G

	K[idofs, idofs] += lc.C'*K[ddofs,idofs] + K[idofs,ddofs]*lc.C + lc.C'*K[ddofs,ddofs]*lc.C
	K[rdofs, idofs] += K[rdofs, ddofs]*lc.C
	K[idofs, rdofs] += lc.C'*K[ddofs, rdofs]

	K[ddofs, rdofs] .= 0.0
	K[ddofs, idofs] .= 0.0
	K[idofs, ddofs] .= 0.0
	K[rdofs, ddofs] .= 0.0
	K[ddofs, ddofs]  = Diagonal(ones(length(ddofs)))
	f[ddofs] .= 0.0
end

function apply_lc!(v, lc)
	v[lc.ddofs] = lc.C*v[lc.idofs] + lc.G
	#@show v[lc.ddofs[1]]
end

