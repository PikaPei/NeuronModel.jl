using CSV
using DataFrames
using Plots


mempots = CSV.read("MemPotALL.dat", delim=' ')
spikes = CSV.read("SpikeALL.dat", delim=' ', header=false)

# neux_spike = spikes[in(["neux"]).(spikes.Column2), :]
# neux_spike = filter(row -> row.Column2 in ["neux"], spikes)  # slower

plot(legend=false)
plot!(mempots[!, 1], mempots[!, 2])
scatter!(spikes[!, 1], fill(35, size(spikes)[1]))