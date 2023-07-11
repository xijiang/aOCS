"""
    function initR(dir)
Initialize R environment for optiSel. These codes only need to be run once.
Call this once before calling `optisel()`.
"""
function initR(dir)
    R"""
        library(optiSel)
        library(data.table)
    """
    ppd = "$dir/ydh"     # phased genotypes for ydh directory
    rbd = "$dir/refs"    # reference breeds directory
    mtd = "$dir/match"   # match directory
    otd = "$dir/other"   # other information directory
    @rput dir ppd rbd mtd otd

    R"""
        animals <- optiSel::read.indiv(file.path(ppd, '0.Chr1.phased'), skip=0, cskip=2)
        bfiles <- paste0(ppd, '/0.Chr', 1:18, '.phased')
        lmp <- data.table::fread(paste0(otd, '/map.txt'))
        fSEG <- segIBD(bfiles, lmp, minSNP = 20, minL = 2.5, keep = animals, skip=0, cskip = 2)

        # calculate native contribution
        Pig <- fread(paste0(otd, '/genotypedindiv.txt'))
        rfiles <- paste0(rbd, '/Others.Chr', 1:18, '.phased')
        mfiles <- paste0(mtd, '/0.Chr', 1:18, '.txt')
        Comp <- optiSel::segBreedComp(mfiles, lmp)
        setnames(Comp, old='native', new='segNC')

        # below needs to be constructed in Julia
        phen <- fread(paste0(otd, '/BV.txt'))
        phen <- merge(phen, Comp[, c("Indiv", "segNC")], on="Indiv")
        phen <- merge(phen, data.frame(Indiv=names(diag(fSEG)), F=(2*fSEG[row(fSEG)==col(fSEG)]-1)), on="Indiv")
        # calculate natKin
        fSEGN <- segIBDatN(list(hap.thisBreed=bfiles, hap.refBreed=rfiles, match=mfiles), Pig, lmp, thisBreed='Y_YDH', minSNP=20, minL=2.5, ubFreq=0.01)
        natKin <- fSEGN$Q1/fSEGN$Q2

        # below needs to be tested
        cont <- data.frame(age=1, male=0.5, female=0.5)
        cand <- candes(phen = phen, fSEG = fSEG, fSEGN = fSEGN, cont = cont)
        Ne <- 100
        L <- 1
        ub.fSEG <- cand$mean$fSEG + (1 - cand$mean$fSEG) / (2 * Ne * L)
        ub.fSEGN <- cand$mean$fSEGN + (1 - cand$mean$fSEGN) / (2 * Ne * L)
        females <- cand$phen$Sex == 'female' & cand$phen$isCandidate
        ub <- setNames(rep(0.00625, sum(females)), cand$phen$Indiv[females])
        con <- list (ub = ub, ub.fSEG = ub.fSEG, ub.fSEGN = ub.fSEGN)
        fit <- opticont('max.segNC', cand, con, solver='cccp', quiet=TRUE)
        Candidate <- fit$parent
        Candidate$n <- noffspring(Candidate, 400, random=TRUE)$nOff
        Mating <- matings(Candidate, Kin=fSEG)
    """
end

"""
    function optisel(dir)
Given the directory `dir` containing the phased genotypes and the map file,
This function finish the R round of optiSel.
"""
function optisel()
    # ToDo: merge genoteyps and origins into one matrix. to feed mating function
    # ToDo: generate a BV file in each generation.
end
