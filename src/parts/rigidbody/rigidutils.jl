
getE(ep::T) where T <: AbstractVector = [-ep[2] ep[1] -ep[4] ep[3]
                                           -ep[3] ep[4] ep[1] -ep[2]
                                           -ep[4] -ep[3] ep[2] ep[1]]

getEbar(ep::T) where T <: AbstractVector = [-ep[2] ep[1] ep[4] -ep[3]
                                                -ep[3] -ep[4] ep[1] ep[2]
                                                -ep[4] ep[3] -ep[2] ep[1]] 
getG(ep::T) where T <: AbstractVector = 2*getE(ep)
getGbar(ep::T) where T <: AbstractVector = 2*getEbar(ep)
getA(ep::T) where T <: AbstractVector = getE(ep)*getEbar(ep)'

#---

getA(t::T) where T <: Real = @SMatrix [cos(t) -sin(t); sin(t) cos(t)]
getA3D(t::T) where T <: Real = @SMatrix [cos(t) -sin(t) T(0); sin(t) cos(t) T(0); T(0) T(0) T(1)]
getdA(t::T) where T <: Real = @SMatrix [-sin(t) -cos(t); cos(t) -sin(t)]
getG(t::T) where T <: Real = @SVector [0, 0, 1]
getGbar(t::T) where T <: Real = @SVector [0, 0, 1]

function remove_unused_cells!(grid, unused_cells::Vector{Int})

    used_nodes = []
    other_cells = setdiff(1:getncells(grid), unused_cells)
    for cellid in other_cells
        nodes = grid.cells[cellid].nodes
        for n in nodes
            push!(used_nodes, n)
        end
    end

    rigid_nodes = []
    for cellid in unused_cells
        nodes = grid.cells[cellid].nodes
        for n in nodes
            push!(rigid_nodes, n)
        end
    end
    rigid_nodes = unique(rigid_nodes)
    allnodes = collect(1:getnnodes(grid))
    unused_nodes = setdiff(allnodes, used_nodes)
    commonnodes = intersect(rigid_nodes, used_nodes)
    #deleteat!(grid.nodes, unused_nodes)
    deleteat!(grid.cells, unused_cells)
    return commonnodes
    #.. delete nodes, delete from cellset etc
end

#=
function rigid_mass_and_intertia(rho::T, grid::Grid{dim}, element::Element, cellset::Vector{Int}) where {T,dim}
    
    #Use the default interpolation
    ip = default_interpolation(celltype(element))
    qr = QuadratureRule{getdim(ip), getrefshape(ip)}(2);
    cv = CellScalarValues(qr,ip)

    mass = 0
    center_of_mass = zero(Vec{dim, Float64})
    inertia = zeros(Float64,3,3)

    #Mass 
    for cell in grid.cells[cellset]
        reinit!(cv, [grid.nodes[id].x for id in cell.nodes])
        for qp in 1:getnquadpoints(cv)
            mass += rho*getdetJdV(cv, qp)
        end
    end

    #Center of mass   
    for cell in grid.cells[cellset]
        reinit!(cv, [grid.nodes[id].x for id in cell.nodes])
        for qp in 1:getnquadpoints(cv)
            xx = zero(Vec{dim, Float64}) 
            for i in 1:getnbasefunctions(cv)
                xx += shape_value(cv, qp, i) * grid.nodes[cell.nodes[i]].x
            end
            center_of_mass += rho*xx*getdetJdV(cv, qp)
        end
    end
    center_of_mass /= mass


    #integrate_mass_andinertia_matrix
    fmi = (rho, x, y) -> rho*(x^2+y^2)
    fpi = (rho,x,y) -> -rho*x*y
    for cell in grid.cells[cellset]#CellsetIterator(dh, element, cellset)
        reinit!(cv, [grid.nodes[id].x for id in cell.nodes])
        for qp in 1:getnquadpoints(cv)
                xx = zero(Vec{dim, Float64})
                
                for iii in 1:getnbasefunctions(cv)
                    xx += shape_value(cv, qp, iii) * grid.nodes[cell.nodes[iii]].x
                end

                xx -= center_of_mass
                xx = Vec{3,T}((xx[1],xx[2], 0.0))

                dΩ = getdetJdV(cv, qp)
                inertia[1,1] += fmi(rho,xx[2],xx[3])*dΩ
                inertia[2,2] += fmi(rho,xx[1],xx[3])*dΩ
                inertia[3,3] += fmi(rho,xx[1],xx[2])*dΩ 

                inertia[1,2] += fpi(rho,xx[1],xx[2])*dΩ 
                inertia[1,3] += fpi(rho,xx[1],xx[3])*dΩ 
                inertia[2,3] += fpi(rho,xx[2],xx[3])*dΩ 

                inertia[2,1] += fpi(rho,xx[2],xx[1])*dΩ 
                inertia[3,1] += fpi(rho,xx[3],xx[1])*dΩ 
                inertia[3,2] += fpi(rho,xx[3],xx[2])*dΩ 
        end

    end

    return mass, center_of_mass, inertia

end
=#
