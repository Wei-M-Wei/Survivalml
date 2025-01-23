SUBROUTINE log_sparse_four_weight (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,pfl1,dfmax,pmax,&
     nlam,flmin,ulam,eps,maxit,intr,nalam,b0,beta,activeGroup,nbeta,alam,&
     npass,jerr,alsparse,lb,ub,weight,intercept_zero)

  USE log_sgl_subfuns_weight
  USE KKT
  IMPLICIT NONE
  ! - - - arg types - - -
  DOUBLE PRECISION, PARAMETER :: big=9.9E30
  DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
  INTEGER, PARAMETER :: mnlam = 6
  INTEGER:: isDifZero
  INTEGER:: mnl
  INTEGER:: bn
  INTEGER:: bs(bn)
  INTEGER:: ix(bn)
  INTEGER:: iy(bn)
  INTEGER:: nobs, nvars, dfmax, pmax, nlam, nalam, npass, jerr, maxit, intr
  INTEGER:: activeGroup(pmax)
  INTEGER:: nbeta(nlam)
  DOUBLE PRECISION :: flmin, eps, alsparse, max_gam, maxDif, al, alf, snorm, d, intercept_zero
  DOUBLE PRECISION, INTENT(in) :: x(nobs,nvars)
  DOUBLE PRECISION, INTENT(in) :: y(nobs)
  DOUBLE PRECISION, INTENT(in) :: weight(nobs)
  DOUBLE PRECISION, INTENT(in) :: pf(bn)
  DOUBLE PRECISION, INTENT(in) :: pfl1(nvars)
  DOUBLE PRECISION :: ulam(nlam)
  DOUBLE PRECISION :: gam(bn)
  DOUBLE PRECISION, INTENT(in) :: lb(bn), ub(bn)
  DOUBLE PRECISION :: beta(nvars,nlam)
  DOUBLE PRECISION :: b0(nlam)
  DOUBLE PRECISION :: alam(nlam)

  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s !need for log_sparse_four_weight
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r ! Residual
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: u

  INTEGER, DIMENSION (:), ALLOCATABLE :: activeGroupIndex
  INTEGER:: g, j, l, ni, me, startix, endix, vl_iter
  DOUBLE PRECISION::t_for_s(bn) ! this is for now just 1/gamma
  ! - - - begin local declarations - - -

  DOUBLE PRECISION:: tlam, lama, lam1ma, al0
  INTEGER:: violation
  INTEGER:: is_in_E_set(bn)
  INTEGER:: is_in_S_set(bn) ! this is for 4-step alg
  DOUBLE PRECISION:: ga(bn)
  DOUBLE PRECISION:: vl(nvars)
  ! - - - allocate variables - - -
  ALLOCATE(b(0:nvars))
  ALLOCATE(oldbeta(0:nvars))
  ALLOCATE(r(1:nobs))
  ALLOCATE(activeGroupIndex(1:bn))
  !    ALLOCATE(al_sparse)
  ! - - - checking pf - ! pf is the relative penalties for each group
  IF(MAXVAL(pf) <= 0.0D0) THEN
     jerr = 10000
     RETURN
  ENDIF
  ! - - - some initial setup - - -
  is_in_E_set = 0
  PRINT *, is_in_E_set
  is_in_S_set = 0
  al = 0.0D0
  mnl = MIN(mnlam, nlam)
  r = 0.0D0
  b = 0.0D0
  oldbeta = 0.0D0
  activeGroup = 0
  activeGroupIndex = 0
  npass = 0
  ni = 0
  alf = 0.0D0
  max_gam = MAXVAL(gam)
  t_for_s = 1 / gam
  IF (intr .ne. 0) THEN
     b(0) = intercept_zero
     r = r + b(0)
  ENDIF
  ! --------- lambda loop ----------------------------
  IF (flmin < 1.0D0) THEN ! THIS is the default...
     flmin = MAX(mfl, flmin) ! just sets a threshold above zero
     alf = flmin ** (1.0D0 / (nlam - 1.0D0))
  ENDIF
  ! PRINT *, alf
  vl = MATMUL(weight*y - 1/(1.0D0+exp(-r)), x)/nobs
  al0 = 0.0D0
  DO g = 1, bn ! For each group...
     ALLOCATE(u(bs(g)))
     u = vl(ix(g):iy(g))
     ga(g) = SQRT(DOT_PRODUCT(u, u))
     DEALLOCATE(u)
  ENDDO
  ! CALL rchkusr()
  DO vl_iter = 1, nvars
     al0 = MAX(al0, ABS(vl(vl_iter))) ! Infty norm of X'y, big overkill for lam_max
  ENDDO
  ! PRINT *, alsparse
  al = al0 ! this value ensures all betas are 0
  l = 0
  tlam = 0.0D0
  DO WHILE (l < nlam)
     ! PRINT *, "Lambda iter =  and al = "
     ! PRINT *, l
     ! PRINT *, al
     ! CALL rchkusr()
     al0 = al
     IF (flmin >= 1.0D0) THEN
        l = l + 1
        al = ulam(l)
     ELSE
        IF (l > 1) THEN
           al = al * alf
           tlam = MAX((2.0*al-al0), 0.0D0)
           l = l + 1
        ELSE IF(l==0) THEN
           al= al * 0.99
           tlam = al
        ENDIF
     ENDIF
     lama = al * alsparse
     lam1ma = al * (1 - alsparse)
     ! This is the start of the algorithm, for a given lambda...
     CALL strong_rule (is_in_S_set, ga, pf, tlam, alsparse)
     ! uses s_set instead of e_set...
     ! --------- outer loop ---------------------------- !
     DO
        ! CALL rchkusr()
        oldbeta(0) = b(0)
        IF (ni > 0) THEN
           DO j = 1, ni
              g = activeGroup(j)
              oldbeta(ix(g):iy(g)) = b(ix(g):iy(g))
           ENDDO
        ENDIF
        ! --inner loop-------------------------------------
        DO
           ! print *, "This is where we enter the inner loop"
           ! CALL rchkusr()
           npass = npass + 1
           maxDif = 0.0D0
           isDifZero = 0 !Boolean to check if b-oldb nonzero. Unnec, in fn.
           DO g = 1, bn
              IF (is_in_E_set(g) == 0) CYCLE
              startix = ix(g)
              endix = iy(g)
              CALL log_update_step_weight(bs(g), startix, endix, b, lama, t_for_s(g), pf(g),&
                   pfl1(startix:endix), lam1ma, x,y, isDifZero, nobs, r, gam(g), maxDif, nvars,&
                   lb(g), ub(g), weight)
              IF (activeGroupIndex(g) == 0 .AND. isDifZero == 1) THEN
                 ni = ni + 1
                 IF (ni > pmax) EXIT
                 activeGroupIndex(g) = ni
                 activeGroup(ni) = g
              ENDIF
           ENDDO
           IF (intr /= 0) THEN
              d = sum(weight*y - 1.0D0/(1.0D0+exp(-r)))
              d = d/nobs
              IF (d /= 0.0D0) THEN
                 b(0) = b(0) + d
                 r= r + d
                 maxDif=MAX(maxDif,d**2)
              ENDIF
           ENDIF
           IF (ni > pmax) EXIT
           IF (maxDif < eps) EXIT
           IF (npass > maxit) THEN
              jerr = -l
              RETURN
           ENDIF
        ENDDO ! End inner loop
        IF (ni > pmax) EXIT
        !--- final check ------------------------
        ! This checks which violate KKT condition
        ! PRINT *, "Here is where the final check starts"
        ! print *, i ! Just to check how many final checks...
        ! i = i+1
        violation = 0
        ! PRINT *, (max_gam * (b - oldbeta) / (1 + ABS(b)))**2
        ! PRINT *, b
        IF (ANY((max_gam * (b - oldbeta) / (1 + ABS(b)))**2 >= eps)) violation = 1
        ! has beta moved globally
        IF (violation == 1) CYCLE
        CALL strong_kkt_check(is_in_E_set, violation, bn, ix, iy, pf, pfl1, lam1ma, bs,&
             lama, ga, is_in_S_set, x, weight*y- 1/(1.0D0+exp(-r)), nobs, nvars, vl) ! Step 3
        IF (violation == 1) CYCLE
        ! Need to compute vl/ga for the ones that aren't already updated,
        ! before log_kkt_check
        ! CALL rchkusr()
        DO g = 1, bn
           IF (is_in_S_set(g) == 0) THEN
              startix = ix(g)
              endix = iy(g)
              ALLOCATE(s(bs(g)))
              s = MATMUL(weight*y- 1.0D0/(1.0D0+exp(-r)), x(:, startix:endix))/nobs
              vl(startix:endix) = s
              CALL softthresh(s, lama*pfl1(startix:endix), bs(g))
              snorm = SQRT(DOT_PRODUCT(s,s))
              ga(g) = snorm
              DEALLOCATE(s)
           ENDIF
        ENDDO
        CALL kkt_check(is_in_E_set, violation, bn, ix, iy, vl, pf, pfl1, lam1ma,&
             bs, lama, ga, nvars) ! Step 4
        IF (violation == 1) CYCLE
        EXIT
     ENDDO ! Ends outer loop
     !---------- final update variable and save results------------
     ! PRINT *, "Here is where the final update starts"
     IF (l == 0) THEN
        IF (MAXVAL(is_in_E_set) == 0) THEN
           CYCLE ! don't save anything, we're still decrementing lambda
        ELSE
           l = 2
           alam(1) = al / MAX(alf, .99D0) ! store previous, larger value
        ENDIF
     ENDIF
     IF(ni > pmax) THEN
        jerr = -10000 - l
        EXIT
     ENDIF
     IF (ni > 0) THEN
        DO j = 1, ni
           g = activeGroup(j)
           beta(ix(g):iy(g),l) = b(ix(g):iy(g))
        ENDDO
     ENDIF
     b0(l) = b(0)
     nbeta(l) = ni
     alam(l) = al
     nalam = l
     IF (l < mnl) CYCLE
     me = 0
     DO j = 1, ni
        g = activeGroup(j)
        IF (ANY(beta(ix(g):iy(g),l) .ne. 0.0D0)) me=me+1
     ENDDO
     IF (me > dfmax) EXIT
  ENDDO ! end lambda loop
  ! print *, is_in_E_set
  DEALLOCATE(b, oldbeta, r, activeGroupIndex)
  RETURN
END SUBROUTINE log_sparse_four_weight

