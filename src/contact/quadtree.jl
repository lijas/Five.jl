

type QuadTree

	maxobjects::Int
	maxlevels::Int

	level::Int
	objects::Vector
	bounds::Rectangle
	nodes::Vector{QuadTree}
    hasnodes::Bool
end

MAXOBJECTS(q::QuadTree) = q.maxobjects
MAXLEVELS(q::QuadTree) = q.maxlevels

function QuadTree(level::Int, bound::Rectangle)
    return QuadTree(6,5,level,Vector(),bound,Vector{QuadTree}(), false)
end

function clear_tree!(quadtree::QuadTree)
    empty!(quadtree.objects);
    if quadtree.hasnodes
        for node in quadtree.nodes
            clear_tree!(node)
        end
    end
    empty!(quadtree.nodes);
end

function split!(quadtree::QuadTree)

    w = quadtree.bounds.w/2;
    h = quadtree.bounds.h/2;
    x = quadtree.bounds.x;
    y = quadtree.bounds.y;

    #nodes[1] = QuadTree(quadtree.level+1, Rectangle(x + w,y,w,h))
    #nodes[2] = QuadTree(quadtree.level+1, Rectangle(x,y,w,h))
    #nodes[3] = QuadTree(quadtree.level+1, Rectangle(x,y+h,w,h))
    #nodes[4] = QuadTree(quadtree.level+1, Rectangle(x + w,y+h,w,h))
    push!(quadtree.nodes, QuadTree(quadtree.level+1, Rectangle(x+w,y,w,h)))
    push!(quadtree.nodes, QuadTree(quadtree.level+1, Rectangle(x,y,w,h)))
    push!(quadtree.nodes, QuadTree(quadtree.level+1, Rectangle(x,y+h,w,h)))
    push!(quadtree.nodes, QuadTree(quadtree.level+1, Rectangle(x+w,y+h,w,h)))
    quadtree.hasnodes = true;
end

function get_quadtree_index(quadtree::QuadTree, object)
    index = -1;

    #Temporarycode
    y = object.y;
    x = object.x;
    w = 0;
    h = 0;
    #Temporarycode

    verticalmidpoint = quadtree.bounds.x + quadtree.bounds.w / 2
    horizontalmidpoint = quadtree.bounds.y + quadtree.bounds.h / 2

    topquadrant = (y< horizontalmidpoint && y+h< horizontalmidpoint)
    bottomquadrant = y > horizontalmidpoint;

    if x < verticalmidpoint && x+w < verticalmidpoint
        if topquadrant
            index = 2;
        elseif bottomquadrant
            index = 3;
        end
    elseif x > verticalmidpoint
        if topquadrant
            index = 1;
        elseif bottomquadrant
            index = 4;
        end
    end

    return index;
end

function insert!(quadtree::QuadTree, object)

    if quadtree.hasnodes
        index = get_quadtree_index(quadtree, object);
        if index != -1
            insert!(quadtree.nodes[index],object)
            return;
        end
    end

    push!(quadtree.objects, object)

    if(length(quadtree.objects) > MAXOBJECTS(quadtree) && quadtree.level < MAXLEVELS(quadtree))
        if !quadtree.hasnodes
            split!(quadtree);
        end

        newobjectlist = [];
        for object in quadtree.objects
            index = get_quadtree_index(quadtree, object)
            if index != -1
                insert!(quadtree.nodes[index], object);
            else
                push!(newobjectlist,object);
            end
        end
        quadtree.objects = newobjectlist;
    end

end

function retrieve(quadtree::QuadTree, returnobjects::Vector, object)
    index = get_quadtree_index(object)
    if index != -1 && quadtree.hasnodes
        retrieve(quadtree.nodes[index], returnobjects, object)
    end

    #add all
    append!(returnobjects, quadtree.objects)

    return returnobjects;
end

function retrieve_bounds(quadtree::QuadTree, returnrecs::Vector)
    
    if quadtree.hasnodes
        for node in quadtree.nodes
            retrieve_bounds(node,returnrecs)
        end
    end

    push!(returnrecs,quadtree.bounds);

    return returnrecs;
end

