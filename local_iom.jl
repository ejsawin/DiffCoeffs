function local_iom(rr,EE,LL)
    v_step=sqrt((v_r(rr,EE,LL)) ^ 2 + (v_t(rr,EE,LL)) ^ 2)

    vp,vp2,vt2=local_coef(rr,EE,LL)

    del_E=0.5*vp2+0.5*vt2+v_step*vp 
    del_E2=(v_step^2)*vp2

    del_L=(LL/v_step)*vp+((rr^2)/(4*LL))*vt2
    del_L2=((LL^2)/(v_step^2))*vp2+0.5*(rr^2-(LL/v_step)^2)*vt2

    del_EL=LL*vp2

    return del_E, del_E2, del_L, del_L2, del_EL
end 