using DelimitedFiles


function output_potential(filename::String, potential)
    potential = round.(potential, digits=3)
    writedlm(filename,  potential, ' ')
end


function output_spike(filename::String, spikes, idx; name=false)
	f = open(filename, "w")
	for spike in spikes
		if name == true
			println(f, spike.second, " ", spike.first, " ", idx[spike.first])
		else
			println(f, spike.second, " ", idx[spike.first])
		end
	end
	close(f)
end
