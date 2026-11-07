# Print the datasets of a checkpoint file with their sizes, without reading array data.
#   julia --project=. regrid/list_checkpoint.jl <checkpoint.jld2>
using JLD2

function walk(file, path)
    object = JLD2.get_dataset(file, path)
    if object.datatype isa Nothing || isempty(path) || object isa JLD2.Group
        return
    end
end

jldopen(ARGS[1], "r") do file
    function visit(group, prefix)
        for key in keys(group)
            path = prefix * key
            child = try
                JLD2.isgroup(file, JLD2.lookup_offset(group, key))
            catch
                false
            end
            if child
                visit(group[key], path * "/")
            else
                dataset = JLD2.get_dataset(group, key)
                println(path, "  ", try JLD2.dims(dataset) catch; try size(dataset.dataspace) catch; "" end end)
            end
        end
    end
    visit(file.root_group, "")
end
