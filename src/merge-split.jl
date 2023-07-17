"""
    function mergegt(dir; nchr = 18, native = "ydh")
This function is to merge the genotypes and breed info on all chromosomes into
one file. Returns a genotype matrix whose element signs indicates the SNP allele
types. The element magnitudes indicates the breed types. This function also
returns two dictionaries, one for converting breed types to numbers, the other
is to convert numbers to breed types.
"""
function mergegt(dir; nchr = 18, native = "ydh")
    pool = joinpath(dir, "pool")
    isdir(pool) || mkdir(pool)

    # Read breed info
    brd = Char[]
    for chr in 1:nchr
        chrfile = joinpath(dir, "match/0.Chr$chr.txt")
        open(chrfile, "r") do io
            readline(io)
            for line in eachline(io)
                for x in split(line)[2:end]
                    push!(brd, x[1])
                end
            end
        end
    end
    dica = Dict{Char, Int8}()
    dica['1'], inc = 1, 2
    for x in unique(brd)
        x == '1' && continue
        dica[x] = inc
        inc += 1
    end
    dicb = Dict{Int8, Char}()
    for (k, v) in dica
        dicb[v] = k
    end
    # Read genotypes
    gt = Int8[]
    mkr, id = String[], nothing
    for chr in 1:nchr
        chrfile = joinpath(dir, native, "0.Chr$chr.phased")
        open(chrfile, "r") do io
            line = readline(io)
            if chr == 1
                id = split(line)[3:2:end]
            end

            for line in eachline(io)
                field = split(line)
                push!(mkr, field[2])
                append!(gt, parse.(Int8, field[3:end]))
            end
        end
    end
    # The mergeing process
    mgt = similar(gt)
    for i in eachindex(gt)
        mgt[i] = gt[i] == 0 ? -dica[brd[i]] : dica[brd[i]]
    end
    # Order the genotypes into a matrix
    # to be finished
    reshape(mgt, 2length(id), :)', id, mkr, dicb
end

"""
    function splitgt(mgt, mkr, dic, lms, grt, dir, native)
Split the merged genotype matrix into two parts, one is the SNP genotypes of 0
and 1. The other is the breed info. Generation `grt` is to be used as the prefix
of the output files. It is also used as the prefix of the ID names in this
generation. The results are saved in `dir`.
"""
function splitgt(mgt, mkr, dic, lms, grt, dir, native)
    imk = 1
    for chr in 1:lms.chr[end]
        open("$dir/$native/$grt.Chr$chr.phased", "w") do io
            print(io, "I IID")
            for id in 1:(size(mgt, 2) ÷ 2)
                print(io, " $grt-$id $grt-$id")
            end
            println(io)
            for _ in 1:lms.nlc[chr]
                print(io, "M ", mkr[imk], ' ')
                println(io, join(Int.(mgt[imk, :] .> 0), ' '))
                imk += 1
            end
        end
    end
    imk = 1
    for chr in 1:lms.chr[end]
        open("$dir/match/$grt.Chr$chr.txt", "w") do io
            print(io, "Name")
            for id in 1:(size(mgt, 2) ÷ 2)
                print(io, " $grt-$id $grt-$id")
            end
            println(io)
            for _ in 1:lms.nlc[chr]
                print(io, mkr[imk])
                for x in mgt[imk, :]
                    print(io, ' ', dic[abs(x)])
                end
                println(io)
                imk += 1
            end
        end
    end
end
