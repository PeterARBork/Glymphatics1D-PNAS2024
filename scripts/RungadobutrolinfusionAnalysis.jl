
using Glymphatics1D

parseandanalyzeGadobutrolinfusion(
    (estimate_D=true, estimate_Dm=false, estimate_v=false),
)

parseandanalyzeGadobutrolinfusion(
    (estimate_D=true, estimate_Dm=true, estimate_v=false),
)

parseandanalyzeGadobutrolinfusion(
    (estimate_D=true, estimate_Dm=false, estimate_v=true),
)

recordings = [
    "210215_mouse03_TACs",
    "210314_mouse02_TACs",
    "210314_mouse03_TACs",
    "210322_mouse04_TACs",
    "210330_mouse05_TACs",
    "210330_mouse06_TACs",
    "210330_mouse07_TACs",
]

indx = 1
raw_datafile = "data/inputs/gadobutrol intrastriatal/$(recordings[indx]).csv"
@time parseandanalyzeGadobutrolinfusion(
    (estimate_D=true, estimate_Dm=true, estimate_v=true, estimate_all_phys=false);
    raw_datafile=raw_datafile,
)

for indx ∈ 1:7
    raw_datafile = "data/inputs/gadobutrol intrastriatal/$(recordings[indx]).csv"
    parseandanalyzeGadobutrolinfusion(
        (estimate_D=true, estimate_Dm=true, estimate_v=true, estimate_all_phys=false);
        raw_datafile=raw_datafile,
    )
end

for indx ∈ 1:7
    raw_datafile = "data/inputs/gadobutrol intrastriatal/$(recordings[indx]).csv"
    @time parseandanalyzeGadobutrolinfusion(
        (estimate_D=true, estimate_Dm=true, estimate_v=true, estimate_all_phys=true);
        raw_datafile=raw_datafile,
    )
end
