using DelimitedFiles
using Plots
using LaTeXStrings

p = [12,13,15,17,19,22,25,27,30,35]

normal = readdlm("rich_mat_normal.txt")
loose = readdlm("rich_mat_loose.txt")
severe = readdlm("rich_mat_severe.txt")
nor_err = readdlm("rich_err_normal.txt")
loo_err = readdlm("rich_err_loose.txt")
sev_err = readdlm("rich_err_severe.txt")

nlp1 = zeros(10,9)
nlp1e = zeros(10,9)
nsp1 = zeros(10,9)
nsp1e = zeros(10,9)
nlp2 = zeros(10,9)
nlp2e = zeros(10,9)
nsp2 = zeros(10,9)
nsp2e = zeros(10,9)
nlm1 = zeros(10,9)
nlm1e = zeros(10,9)
nsm1 = zeros(10,9)
nsm1e = zeros(10,9)
nlm2 = zeros(10,9)
nlm2e = zeros(10,9)
nsm2 = zeros(10,9)
nsm2e = zeros(10,9)
nnp1 = zeros(10,9)
nnp1e = zeros(10,9)
nnn1 = zeros(10,9)
nnn1e = zeros(10,9)
nK1p = zeros(10)
nK2p = zeros(10)
nK3p = zeros(10)
eK1p = zeros(10)
eK2p = zeros(10)
eK3p = zeros(10)
nK1m = zeros(10)
nK2m = zeros(10)
nK3m = zeros(10)
eK1m = zeros(10)
eK2m = zeros(10)
eK3m = zeros(10)
outputp = zeros(10,6)
outputm = zeros(10,6)

for i in 1:10
    for j in 1:9
        nnp1[i,div(j-1,3)+mod(j-1,3)*3+1] = normal[2*(i+4)-1,2+div(j+2,3)*3+j]
        nnp1e[i,div(j-1,3)+mod(j-1,3)*3+1] = nor_err[2*(i+4)-1,2+div(j+2,3)*3+j]
        nnn1[i,div(j-1,3)+mod(j-1,3)*3+1] = normal[2*(i+4)-1,2+div(j-1,3)*3+j]
        nnn1e[i,div(j-1,3)+mod(j-1,3)*3+1] = nor_err[2*(i+4)-1,2+div(j-1,3)*3+j]
        nlp1[i,div(j-1,3)+mod(j-1,3)*3+1] = normal[2*(i+4)-1,2+div(j+2,3)*3+j]-loose[2*(i+4)-1,2+div(j+2,3)*3+j]
        nlp1e[i,div(j-1,3)+mod(j-1,3)*3+1] = sqrt(nor_err[2*(i+4)-1,2+div(j+2,3)*3+j]^2+loo_err[2*(i+4)-1,2+div(j+2,3)*3+j]^2)
        # nlp2[i,div(j-1,3)+mod(j-1,3)*3+1] = normal[2*(i+4),2+div(j+2,3)*3+j]-loose[2*(i+4),2+div(j+2,3)*3+j]
        # nlp2e[i,div(j-1,3)+mod(j-1,3)*3+1] = sqrt(nor_err[2*(i+4),2+div(j+2,3)*3+j]^2+loo_err[2*(i+4),2+div(j+2,3)*3+j]^2)
        nsp1[i,div(j-1,3)+mod(j-1,3)*3+1] = normal[2*(i+4)-1,2+div(j+2,3)*3+j]-severe[2*(i+4)-1,2+div(j+2,3)*3+j]
        nsp1e[i,div(j-1,3)+mod(j-1,3)*3+1] = sqrt(nor_err[2*(i+4)-1,2+div(j+2,3)*3+j]^2+sev_err[2*(i+4)-1,2+div(j+2,3)*3+j]^2)
        # nsp2[i,div(j-1,3)+mod(j-1,3)*3+1] = normal[2*(i+4),2+div(j+2,3)*3+j]-severe[2*(i+4),2+div(j+2,3)*3+j]
        # nsp2e[i,div(j-1,3)+mod(j-1,3)*3+1] = sqrt(nor_err[2*(i+4),2+div(j+2,3)*3+j]^2+sev_err[2*(i+4),2+div(j+2,3)*3+j]^2)
        nlm1[i,div(j-1,3)+mod(j-1,3)*3+1] = normal[2*(i+4)-1,2+div(j-1,3)*3+j]-loose[2*(i+4)-1,2+div(j-1,3)*3+j]
        nlm1e[i,div(j-1,3)+mod(j-1,3)*3+1] = sqrt(nor_err[2*(i+4)-1,2+div(j-1,3)*3+j]^2+loo_err[2*(i+4)-1,2+div(j-1,3)*3+j]^2)
        # nlm2[i,div(j-1,3)+mod(j-1,3)*3+1] = normal[2*(i+4),2+div(j-1,3)*3+j]-loose[2*(i+4),2+div(j-1,3)*3+j]
        # nlm2e[i,div(j-1,3)+mod(j-1,3)*3+1] = sqrt(nor_err[2*(i+4),2+div(j-1,3)*3+j]^2+loo_err[2*(i+4),2+div(j-1,3)*3+j]^2)
        nsm1[i,div(j-1,3)+mod(j-1,3)*3+1] = normal[2*(i+4)-1,2+div(j-1,3)*3+j]-severe[2*(i+4)-1,2+div(j-1,3)*3+j]
        nsm1e[i,div(j-1,3)+mod(j-1,3)*3+1] = sqrt(nor_err[2*(i+4)-1,2+div(j-1,3)*3+j]^2+sev_err[2*(i+4)-1,2+div(j-1,3)*3+j]^2)
        # nsm2[i,div(j-1,3)+mod(j-1,3)*3+1] = normal[2*(i+4),2+div(j-1,3)*3+j]-severe[2*(i+4),2+div(j-1,3)*3+j]
        # nsm2e[i,div(j-1,3)+mod(j-1,3)*3+1] = sqrt(nor_err[2*(i+4),2+div(j-1,3)*3+j]^2+sev_err[2*(i+4),2+div(j-1,3)*3+j]^2)
    end
    nK1p[i] = nnp1[i,2]
    nK2p[i] = nnp1[i,5]
    nK3p[i] = nnp1[i,8]
    eK1p[i] = nnp1e[i,2]
    eK2p[i] = nnp1e[i,5]
    eK3p[i] = nnp1e[i,8]
    nK1m[i] = nnn1[i,2]
    nK2m[i] = nnn1[i,5]
    nK3m[i] = nnn1[i,8]
    eK1m[i] = nnn1e[i,2]
    eK2m[i] = nnn1e[i,5]
    eK3m[i] = nnn1e[i,8]
    outputp[i,1] = nK1p[i]
    outputp[i,2] = eK1p[i]
    outputp[i,3] = nK2p[i]
    outputp[i,4] = eK2p[i]
    outputp[i,5] = nK3p[i]
    outputp[i,6] = eK3p[i]
    outputm[i,1] = nK1m[i]
    outputm[i,2] = eK1m[i]
    outputm[i,3] = nK2m[i]
    outputm[i,4] = eK2m[i]
    outputm[i,5] = nK3m[i]
    outputm[i,6] = eK3m[i]
end

plot(p,nK1p,seriestype=:scatter, ylims = (-0.1,1.1), marker = (:circle), markercolor = :blue, yerror=eK1p, label=L"\pi^+ \rightarrow K^+", guidefont=font(15), tickfont=font(12),
        xlabel = L"p_{h}\;(GeV/c)", ylabel = L"P(h^{+} \rightarrow K^{+})", legend=:inside, legendfontsize=12)
plot!(p,nK2p,seriestype=:scatter, ylims = (-0.1,1.1), marker = (:diamond), markercolor = :green, yerror=eK2p, label=L"K^+ \rightarrow K^+", guidefont=font(15), tickfont=font(12),
        xlabel = L"p_{h}\;(GeV/c)", ylabel = L"P(h^{+} \rightarrow K^{+})", legend=:inside, legendfontsize=12)
plot!(p,nK3p,seriestype=:scatter, ylims = (-0.1,1.1), marker = (:square), markercolor = :red, yerror=eK3p, label=L"p \rightarrow K^+", guidefont=font(15), tickfont=font(12),
        xlabel = L"p_{h}\;(GeV/c)", ylabel = L"P(h^{+} \rightarrow K^{+})", legend=:inside, legendfontsize=12)
savefig("RICHKp.png")

plot(p,nK1m,seriestype=:scatter, ylims = (-0.1,1.1), marker = (:circle), markercolor = :blue, yerror=eK1m, label=L"\pi^- \rightarrow K^-", guidefont=font(15), tickfont=font(12),
        xlabel = L"p_{h}\;(GeV/c)", ylabel = L"P(h^{-} \rightarrow K^{-})", legend=:inside, legendfontsize=12)
plot!(p,nK2m,seriestype=:scatter, ylims = (-0.1,1.1), marker = (:diamond), markercolor = :green, yerror=eK2m, label=L"K^- \rightarrow K^-", guidefont=font(15), tickfont=font(12),
        xlabel = L"p_{h}\;(GeV/c)", ylabel = L"P(h^{-} \rightarrow K^{-})", legend=:inside, legendfontsize=12)
plot!(p,nK3m,seriestype=:scatter, ylims = (-0.1,1.1), marker = (:square), markercolor = :red, yerror=eK3m, label=L"\bar{p} \rightarrow K^-", guidefont=font(15), tickfont=font(12),
        xlabel = L"p_{h}\;(GeV/c)", ylabel = L"P(h^{-} \rightarrow K^{-})", legend=:inside, legendfontsize=12)
savefig("RICHKm.png")

plot(p,nnp1,layout=(3,3),seriestype=:scatter, size = (1500,1500), dpi=300, ylims = (0,1), marker = (:circle), markercolor = :blue, yerror=nnp1e, legend=false, guidefont=font(15),
        xlabel = L"p_{h}\;(GeV/c)", ylabel = [L"\Delta P(\pi^{+} \rightarrow \pi^{+})" L"\Delta P(\pi^{+}\rightarrow K^{+})" L"\Delta P(\pi^{+} \rightarrow p)" L"\Delta P(K^{+} \rightarrow \pi^{+})" L"\Delta P(K^{+}\rightarrow K^{+})" L"\Delta P(K^{+} \rightarrow p)" L"\Delta P(p \rightarrow \pi^{+})" L"\Delta P(p \rightarrow K^{+})" L"\Delta P(p \rightarrow p)"])
savefig("RICHPlus.png")

plot(p,nnn1,layout=(3,3),seriestype=:scatter, size = (1500,1500), dpi=300, ylims = (0,1), marker = (:circle), markercolor = :blue, yerror=nnn1e, legend=false, guidefont=font(15),
        xlabel = L"p_{h}\;(GeV/c)", ylabel = [L"\Delta P(\pi^{+} \rightarrow \pi^{+})" L"\Delta P(\pi^{+}\rightarrow K^{+})" L"\Delta P(\pi^{+} \rightarrow p)" L"\Delta P(K^{+} \rightarrow \pi^{+})" L"\Delta P(K^{+}\rightarrow K^{+})" L"\Delta P(K^{+} \rightarrow p)" L"\Delta P(p \rightarrow \pi^{+})" L"\Delta P(p \rightarrow K^{+})" L"\Delta P(p \rightarrow p)"])
savefig("RICHMinus.png")

plot(p,nlp1,layout=(3,3),seriestype=:scatter, size = (1500,1500), dpi=300, ylims = (-0.1,0.1),legend=:bottom, marker = (:circle), markercolor = :green, yerror=nlp1e, guidefont=font(15), label="Loose",
        xlabel = L"p_{h}\;(GeV/c)", ylabel = [L"\Delta P(\pi^{+} \rightarrow \pi^{+})" L"\Delta P(\pi^{+}\rightarrow K^{+})" L"\Delta P(\pi^{+} \rightarrow p)" L"\Delta P(K^{+} \rightarrow \pi^{+})" L"\Delta P(K^{+}\rightarrow K^{+})" L"\Delta P(K^{+} \rightarrow p)" L"\Delta P(p \rightarrow \pi^{+})" L"\Delta P(p \rightarrow K^{+})" L"\Delta P(p \rightarrow p)"])
# plot!(p,nlp2,layout=(3,3),seriestype=:scatter, size = (1500,1500), dpi=300, ylims = (-0.1,0.1),legend=false, marker = (:circle), markercolor = :green, yerror=nlp2e,guidefont=font(15),
        # xlabel = L"p_{h}\;(GeV/c)", ylabel = [L"\Delta P(\pi^{+} \rightarrow \pi^{+})" L"\Delta P(\pi^{+}\rightarrow K^{+})" L"\Delta P(\pi^{+} \rightarrow p)" L"\Delta P(K^{+} \rightarrow \pi^{+})" L"\Delta P(K^{+}\rightarrow K^{+})" L"\Delta P(K^{+} \rightarrow p)" L"\Delta P(p \rightarrow \pi^{+})" L"\Delta P(p \rightarrow K^{+})" L"\Delta P(p \rightarrow p)"])
plot!(p,nsp1,layout=(3,3),seriestype=:scatter, size = (1500,1500), dpi=300, ylims = (-0.1,0.1), legend=:bottom, marker = (:square), markercolor = :blue, yerror=nsp1e,guidefont=font(15), label="Severe",
        xlabel = L"p_{h}\;(GeV/c)", ylabel = [L"\Delta P(\pi^{+} \rightarrow \pi^{+})" L"\Delta P(\pi^{+}\rightarrow K^{+})" L"\Delta P(\pi^{+} \rightarrow p)" L"\Delta P(K^{+} \rightarrow \pi^{+})" L"\Delta P(K^{+}\rightarrow K^{+})" L"\Delta P(K^{+} \rightarrow p)" L"\Delta P(p \rightarrow \pi^{+})" L"\Delta P(p \rightarrow K^{+})" L"\Delta P(p \rightarrow p)"])
# plot!(p,nsp2,layout=(3,3),seriestype=:scatter, size = (1500,1500), dpi=300, ylims = (-0.1,0.1), legend=false, marker = (:square), markercolor = :green, yerror=nsp2e,guidefont=font(15),
        # xlabel = L"p_{h}\;(GeV/c)", ylabel = [L"\Delta P(\pi^{+} \rightarrow \pi^{+})" L"\Delta P(\pi^{+}\rightarrow K^{+})" L"\Delta P(\pi^{+} \rightarrow p)" L"\Delta P(K^{+} \rightarrow \pi^{+})" L"\Delta P(K^{+}\rightarrow K^{+})" L"\Delta P(K^{+} \rightarrow p)" L"\Delta P(p \rightarrow \pi^{+})" L"\Delta P(p \rightarrow K^{+})" L"\Delta P(p \rightarrow p)"])
savefig("SysPlus.png")

plot(p,nlm1,layout=(3,3),seriestype=:scatter, size = (1500,1500), dpi=300, ylims = (-0.1,0.1),legend=:bottom, marker = (:circle), markercolor = :green, yerror=nlm1e,guidefont=font(15), label="Loose",
        xlabel = L"p_{h}\;(GeV/c)", ylabel = [L"\Delta P(\pi^{-} \rightarrow \pi^{-})" L"\Delta P(\pi^{-}\rightarrow K^{-})" L"\Delta P(\pi^{-} \rightarrow \bar{p})" L"\Delta P(K^{-} \rightarrow \pi^{-})" L"\Delta P(K^{-}\rightarrow K^{-})" L"\Delta P(K^{-} \rightarrow \bar{p})" L"\Delta P(\bar{p} \rightarrow \pi^{-})" L"\Delta P(\bar{p} \rightarrow K^{-})" L"\Delta P(\bar{p} \rightarrow \bar{p})"])
# plot!(p,nlm2,layout=(3,3),seriestype=:scatter, size = (1500,1500), dpi=300, ylims = (-0.1,0.1),legend=false, marker = (:circle), markercolor = :green, yerror=nlm2e,guidefont=font(15),
        # xlabel = L"p_{h}\;(GeV/c)", ylabel = [L"\Delta P(\pi^{-} \rightarrow \pi^{-})" L"\Delta P(\pi^{-}\rightarrow K^{-})" L"\Delta P(\pi^{-} \rightarrow \bar{p})" L"\Delta P(K^{-} \rightarrow \pi^{-})" L"\Delta P(K^{-}\rightarrow K^{-})" L"\Delta P(K^{-} \rightarrow \bar{p})" L"\Delta P(\bar{p} \rightarrow \pi^{-})" L"\Delta P(\bar{p} \rightarrow K^{-})" L"\Delta P(\bar{p} \rightarrow \bar{p})"])
plot!(p,nsm1,layout=(3,3),seriestype=:scatter, size = (1500,1500), dpi=300, ylims = (-0.1,0.1), legend=:bottom, marker = (:square), markercolor = :blue, yerror=nsm1e,guidefont=font(15), label="Severe",
        xlabel = L"p_{h}\;(GeV/c)", ylabel = [L"\Delta P(\pi^{-} \rightarrow \pi^{-})" L"\Delta P(\pi^{-}\rightarrow K^{-})" L"\Delta P(\pi^{-} \rightarrow \bar{p})" L"\Delta P(K^{-} \rightarrow \pi^{-})" L"\Delta P(K^{-}\rightarrow K^{-})" L"\Delta P(K^{-} \rightarrow \bar{p})" L"\Delta P(\bar{p} \rightarrow \pi^{-})" L"\Delta P(\bar{p} \rightarrow K^{-})" L"\Delta P(\bar{p} \rightarrow \bar{p})"])
# plot!(p,nsm2,layout=(3,3),seriestype=:scatter, size = (1500,1500), dpi=300, ylims = (-0.1,0.1), legend=false, marker = (:square), markercolor = :green, yerror=nsm2e,guidefont=font(15),
        # xlabel = L"p_{h}\;(GeV/c)", ylabel = [L"\Delta P(\pi^{-} \rightarrow \pi^{-})" L"\Delta P(\pi^{-}\rightarrow K^{-})" L"\Delta P(\pi^{-} \rightarrow \bar{p})" L"\Delta P(K^{-} \rightarrow \pi^{-})" L"\Delta P(K^{-}\rightarrow K^{-})" L"\Delta P(K^{-} \rightarrow \bar{p})" L"\Delta P(\bar{p} \rightarrow \pi^{-})" L"\Delta P(\bar{p} \rightarrow K^{-})" L"\Delta P(\bar{p} \rightarrow \bar{p})"])
savefig("SysMinus.png")

writedlm("KplusEff.txt", outputp)
writedlm("KminusEff.txt", outputm)
