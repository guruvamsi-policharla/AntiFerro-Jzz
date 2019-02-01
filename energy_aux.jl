function total_mag(lat)
    #returns total magnetisation vector
    return sum(lat)
end

function test_flip(x, y, J, lat, T)
""" Checks whether energy allows for a flip or not """
    #a = sample_uni()
    a = sample_gauss(lat[x,y])
    de = -energy_pos(x,y,J,lat) + energy_pos(x,y,J,lat,a);
    
    if(de<0)
        lat[x,y] = a
        return true
    elseif(rand() < exp(-de/T))
        lat[x,y] = a
        return true
    else
        return false
    end
end

function myplus(a,b,N)
    c = mod(a+b,N)
    if c == 0
        c = N
    end
    return c
end


function energy_pos(x, y, J, lat, a = [0,0,0])
    N = size(lat,1)

    left = lat[myplus(x,-1,N),y]
    right = lat[myplus(x,1,N),y]
    down = lat[x,myplus(y,1,N)]
    up = lat[x,myplus(y,-1,N)]

    if(a == [0,0,0])
        energy = 1*dot(lat[x,y],left + right + up + down) + (J-1)*(lat[x,y][3])*(left[3] + up[3] + right[3] + down[3]);
        return energy
    else
        energy = 1*dot(a,left + right + up + down) + (J-1)*a[3]*(left[3] + up[3] + right[3] + down[3]);
        return energy
    end
end

function total_energy(J,lat)
    e = 0.0
    for i = 1:size(lat,1)
        for j = 1:size(lat,2)
            e = e + energy_pos(i,j,J,lat)
        end
    end
    return e/2;
end
