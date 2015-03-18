# some utilities for parallel processing

function get_chunks(n::Int64)
    chunks_i = round(Int,linspace(1,n,1+length(workers())))
    chunks_f = chunks_i[2:end]
    chunks_f[1:end-1] -= 1
    chunks_i = chunks_i[1:end-1]
    [chunks_i[i]:chunks_f[i] for i in 1:length(chunks_i)]
end
