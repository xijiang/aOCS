"""
    function xps(dir)

## Description

This function do the following jobs:
- Initialize the `R` environment
- In each repeat:
  - find the optimum mating pairs of the current generation with `optiSel`
  - merge genotypes and breed info into one file
  - reprodude the next generation
  - split genotypes and breed info for `optiSel` to work on them
- Summarize the results

## Notes

- All the data are in `dir = "dat"`
- Specify the native breed folder when necessary, default `ydh`.
  - The file names should be `0.Chr1.phased`, `0.Chr2.phased`, etc.
  - Where `0` is the generation number.
- The phased references should be in `refs`.
  - The file names should be `Others.Chr1.phased`, `Others.Chr2.phased`, etc.
  - Above is for the strange requirement of `optiSel`.
- Linkage map, phenotypes, and breeding values are in `other` folder.
- The match results for breed info in the genotypes are in `match`.
  - This can be generated with some other `optiSel` functions
  - Here it is once for all.
"""
function xps(ngrt = 5, noff = 400; dir = "dat", native="ydh")
    mgt, animals, mkr, dic= mergegt(dir, native = native)
    nmk = length(mkr)
    initR(dir, noff)
    @rget lmp
    lms = sumMap(DataFrame(chr=Int8.(lmp.Chr), pos = Int.(floor.(lmp.Mb .* 1e6))))
    for igrt in 0:ngrt
      @info "Generation $igrt"
      replaceid(dir, igrt, animals)
      matings = optipm(igrt)
      pm = pkped(matings, animals)
      nid = size(pm, 1)
      og = zeros(Int8, nmk, 2nid)
      drop(mgt, og, pm, lms)
      mgt = og
      og = nothing
      splitgt(mgt, mkr, dic, lms, igrt+1, dir, native)
      initpt(igrt+1, nid)
      @rget animals
    end
    summarize()
end
