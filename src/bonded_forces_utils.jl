

"""
Return the bead id of the previous bead on a chain

chainbounds is a tuple of at least 2 Ints, 
    chainbounds[1] is 1
    chainbounds[end] is the total number of beads + 1
    the other element in chainbounds are the bead ids of the start of a new chain
"""
function prev_beadid(beadid, chainbounds)
    prev= beadid - 1
    #wrap around for each chain, usually there are only 1 or 2 chains.
    for i in 1:(length(chainbounds)-1)
        if beadid==chainbounds[i] #wrap around to end of chain
            prev= chainbounds[i+1]-1
        end
    end
    prev
end

"""
Return the bead id of the next bead on a chain

chainbounds is a tuple of at least 2 Ints, 
    chainbounds[1] is 1
    chainbounds[end] is the total number of beads + 1
    the other element in chainbounds are the bead ids of the start of a new chain
"""
function next_beadid(beadid, chainbounds)
    next= beadid + 1
    #wrap around for each chain, usually there are only 1 or 2 chains.
    for i in 2:length(chainbounds)
        if next==chainbounds[i] #wrap around to start of chain
            next= chainbounds[i-1]
        end
    end
    next
end
