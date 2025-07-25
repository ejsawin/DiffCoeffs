function average_coef_plot(E_fixed,lower,upper,steps)
    ll_range=exp.(range(log(lower),log(upper),length=steps))

    de=abs.(average_coef.(E_fixed,ll_range,1))
    de2=abs.(average_coef.(E_fixed,ll_range,2))
    dl=abs.(average_coef.(E_fixed,ll_range,3))
    dl2=abs.(average_coef.(E_fixed,ll_range,4))
    del=abs.(average_coef.(E_fixed,ll_range,5))

    graph1=plot(ll_range,[de de2 dl dl2 del],label=["DE" "DE2" "DL" "DL2" "DEL"],title="Omega Cen Diff Coefs",xscale=:log10,yscale=:log10,xlabel="log_10(L)",ylabel="log_10")
    display(graph1)
    savefig("coef_plot.png")
    readline()
end