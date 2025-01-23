MODULE KKT

   USE spmatmul
   IMPLICIT NONE

   CONTAINS

   SUBROUTINE strong_rule (is_in_E_set, ga, pf, tlam, alsparse)
      IMPLICIT NONE
      INTEGER :: g, k
      INTEGER, DIMENSION (:), INTENT(inout) :: is_in_E_set
      DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: ga
      DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: pf
      DOUBLE PRECISION, INTENT(in) :: tlam, alsparse
      DOUBLE PRECISION :: z
      k = SIZE(is_in_E_set)
      z = tlam * (1 - alsparse)

      DO g = 1, k
         IF (is_in_E_set(g) == 1) CYCLE
         IF (ga(g) > pf(g) * z) is_in_E_set(g) = 1
      ENDDO
      RETURN
   END SUBROUTINE strong_rule

   SUBROUTINE kkt_check(is_in_E_set, violation, bn, ix, iy, vl, pf,pfl1,&
        lam1ma, bs, lama, ga, nvars)
      IMPLICIT NONE
      INTEGER :: g, startix, endix, nvars
      INTEGER, INTENT(in) :: bn
      INTEGER, INTENT(in) :: bs(bn)
      INTEGER, INTENT(in) :: ix(bn), iy(bn)
      INTEGER, DIMENSION(:), INTENT(inout) :: is_in_E_set
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: ga
      DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: vl
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s
      DOUBLE PRECISION :: snorm
      DOUBLE PRECISION, INTENT(in) :: pf(bn)
      DOUBLE PRECISION, INTENT(in) :: pfl1(nvars)
      INTEGER, INTENT(inout) :: violation
      DOUBLE PRECISION, INTENT(in) :: lam1ma, lama

      DO g = 1, bn
         IF (is_in_E_set(g) == 1) CYCLE
         startix = ix(g)
         endix = iy(g)
         ALLOCATE(s(bs(g)))
         s = vl(startix:endix)
         CALL softthresh(s, lama * pfl1(startix:endix), bs(g))
         snorm = SQRT(DOT_PRODUCT(s,s))
         ga(g) = snorm
         IF(ga(g) > pf(g) * lam1ma) THEN
            is_in_E_set(g) = 1
            violation = 1
         ENDIF
         DEALLOCATE(s)
      ENDDO
      RETURN
   END SUBROUTINE kkt_check


   SUBROUTINE strong_kkt_check(is_in_E_set,violation,bn,ix,iy,pf,pfl1,lam1ma,bs,&
         lama,ga,is_in_S_set,x,r,nobs,nvars,vl)
      IMPLICIT NONE
      INTEGER, INTENT(in)::nobs
      INTEGER, INTENT(in)::nvars
      DOUBLE PRECISION,INTENT(in):: x(nobs, nvars)
      DOUBLE PRECISION, INTENT(in):: r(nobs)
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: vl
      INTEGER :: g, startix, endix
      INTEGER, INTENT(in) :: bn
      INTEGER, INTENT(in) ::bs(bn)
      INTEGER, INTENT(in) :: ix(bn), iy(bn)
      INTEGER, DIMENSION(:), INTENT(inout) :: is_in_E_set
      INTEGER, DIMENSION(:), INTENT(in) :: is_in_S_set
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: ga
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s
      DOUBLE PRECISION :: snorm
      DOUBLE PRECISION, INTENT(in) :: pf(bn)
      DOUBLE PRECISION, INTENT(in) :: pfl1(nvars)
      INTEGER, INTENT(inout) :: violation
      DOUBLE PRECISION, INTENT(in) :: lam1ma, lama

      violation = 0
      DO g = 1, bn
         IF (is_in_S_set(g) == 1) THEN
            startix = ix(g)
            endix = iy(g)
            ALLOCATE(s(bs(g)))
            s = MATMUL(r, x(:,startix:endix)) / nobs
            vl(startix:endix) = s
            CALL softthresh(s, lama * pfl1(startix:endix), bs(g))
            snorm = SQRT(dot_PRODUCT(s,s))
            ga(g) = snorm
            DEALLOCATE(s)
            IF (is_in_E_set(g) == 1) CYCLE
            IF (ga(g) > pf(g) * lam1ma) THEN
               is_in_E_set(g) = 1
               violation = 1
            ENDIF
         ENDIF
      ENDDO
      RETURN
   END SUBROUTINE strong_kkt_check


END MODULE KKT


