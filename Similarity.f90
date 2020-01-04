!! Coded by zxli 2020.01.01
!! Contact Info : flyaway333f@163.com
!! Reference :
!! https://wenku.baidu.com/view/3da98708876fb84ae45c3b3567ec102de2bddfaa.html
!! https://www.cnblogs.com/liujinhong/p/6001997.html
!! http://www.cppblog.com/unixfy/archive/2012/02/13/165487.html
module Similarity
    implicit none
contains

!! AVG 平均值
!! Variance 方差
!! Standard_Deviation 标准差
!! Covariance 协方差
!! Covariance_Matrix 协方差矩阵
!! Inverse 逆矩阵

!! Euclidean_Distance 欧几里得距离
!! Manhattan_Distance 曼哈顿距离
!! Chebyshev_Distance 切比雪夫距离
!! Minkowski_Distance 闵可夫斯基距离
!! Standardized_Euclidean_Distance 标准化欧几里得距离
!! Mahalanobis_Distance 马氏距离
!! Cosine 夹角余弦
!! Hamming_Distance 汉明距离
!! Jaccard_Similarity_Coefficient 杰卡德相似系数(未实现)
!! Correlation_Coefficient 相关系数
!! Tanimoto_Coefficient Tanimoto系数
!! Lance_Williams_Distance 兰氏距离
!! Dice_Coefficient Dice系数

pure function AVG(n,v)
    implicit none
    real(kind=8)            :: AVG
    integer,intent(in)      :: n
    real(kind=8),intent(in) :: v(n)
    AVG = SUM(v,n)/n
end function AVG

pure function Variance(n,v)
    implicit none
    real(kind=8)            :: Variance
    integer,intent(in)      :: n
    real(kind=8),intent(in) :: v(n)
    real(kind=8) :: mean
    real(kind=8) :: d(n)
    mean = AVG(n,v)
    d = v - mean
    Variance = DOT_PRODUCT(d,d)/(n-1)
end function Variance

pure function Standard_Deviation(n,v)
    implicit none
    real(kind=8)            :: Standard_Deviation
    integer,intent(in)      :: n
    real(kind=8),intent(in) :: v(n)
    Standard_Deviation = DSQRT(Variance(n,v))
end function Standard_Deviation

pure function Covariance(n,u,v)
    implicit none
    real(kind=8)            :: Covariance
    integer,intent(in)      :: n
    real(kind=8),intent(in) :: u(n), v(n)
    real(kind=8) :: mean_u, mean_v
    real(kind=8) :: d_u(n), d_v(n)
    mean_u = AVG(n,u)
    mean_v = AVG(n,v)
    d_u = u - mean_u
    d_v = v - mean_v
    Covariance = DOT_PRODUCT(d_u,d_v)/(n-1)
end function Covariance

subroutine Covariance_Matrix(m,n,X,S)
    implicit none
    !! there are n data with dimension equal m
    integer,intent(in)      :: m,n
    real(kind=8),intent(in) :: X(m,n)
    real(kind=8) :: S(m,m)
    integer      :: i, j
    Forall (i=1:m)
        S(i,i)=Variance(n,X(i,:))
    End Forall
    Forall (i=1:m,j=1:m,i<j)
        S(i,j) = Covariance(n,X(i,:),X(j,:))
    End Forall
    Forall (i=1:m,j=1:m,i>j)
        S(i,j) = S(j,i)
    End Forall
end subroutine Covariance_Matrix

subroutine Inverse(n,M,iM)
    !! Learn from
    !! http://www.yunsuan.info/cgi-bin/matrix_inverse.py
    !! https://blog.csdn.net/lcj_cjfykx/article/details/44803477
    implicit none
    integer,intent(in)       :: n
    real(kind=8),intent(in)  :: M(n,n)
    real(kind=8),intent(out) :: iM(n,n)
    integer      :: ipiv(N)
    integer      :: info1, info2
    integer      :: i
    call DGETRF(n,n,M,N,ipiv,info1)
    if (info1/=0) write(*,"(A)") "DGETRF failed"
    iM = 0.0d0
    Forall (i=1:n)
        iM(i,i) = 1.0d0
    End Forall
    call DGETRS('N',n,n,M,n,ipiv,iM,n,info2)
    if (info2/=0) write(*,"(A)") "DGETRS failed"
end subroutine Inverse

subroutine Euclidean_Distance(n,u,v,dist)
    implicit none
    integer,intent(in)      :: n
    real(kind=8),intent(in) :: u(n), v(n)
    real(kind=8),intent(out):: dist
    real(kind=8)            :: diff(n)
    diff = u - v
    dist = DSQRT(DOT_PRODUCT(diff,diff))
end subroutine Euclidean_Distance

subroutine Manhattan_Distance(n,u,v,dist)
    implicit none
    integer,intent(in)      :: n
    real(kind=8),intent(in) :: u(n), v(n)
    real(kind=8),intent(out):: dist
    real(kind=8)            :: diff(n)
    diff = u - v
    diff = DABS(diff)
    dist = SUM(diff,n)
end subroutine Manhattan_Distance

subroutine Chebyshev_Distance(n,u,v,dist)
    implicit none
    integer,intent(in)      :: n
    real(kind=8),intent(in) :: u(n), v(n)
    real(kind=8),intent(out):: dist
    real(kind=8)            :: diff(n)
    diff = u - v
    diff = DABS(diff)
    dist = MAXVAL(diff,n)
end subroutine Chebyshev_Distance

subroutine Minkowski_Distance(n,u,v,dist,p)
    implicit none
    integer,intent(in)      :: n
    real(kind=8),intent(in) :: u(n), v(n)
    real(kind=8),intent(out):: dist
    integer,intent(in)      :: p
    real(kind=8)            :: diff(n)
    diff = u - v
    diff = DABS(diff)
    diff = diff**p
    dist = SUM(diff,n)
    dist = dist**(1.0/p)
end subroutine Minkowski_Distance

subroutine Standardized_Euclidean_Distance(n,u,v,dist,s)
    implicit none
    integer,intent(in)      :: n
    real(kind=8),intent(in) :: u(n), v(n)
    real(kind=8),intent(out):: dist
    real(kind=8),intent(in) :: s(n)
    real(kind=8)            :: diff(n)
    diff = u - v
    diff = diff/s
    dist = DOT_PRODUCT(diff,diff)
    dist = DSQRT(dist)
end subroutine Standardized_Euclidean_Distance

subroutine Mahalanobis_Distance(m,n,X,D)
    !! Learn from
    !! https://blog.csdn.net/lcj_cjfykx/article/details/44803477
    !! https://cloud.tencent.com/developer/article/1447079
    implicit none
    !! there are n data with dimension equal m
    integer,intent(in)      :: m,n
    real(kind=8),intent(in) :: X(m,n)
    real(kind=8),intent(out):: D(n,n)
    real(kind=8) :: S(m,m), iS(m,m)
    real(kind=8) :: Delta(1,m,n,n)
    real(kind=8) :: Temp (1,m,n,n)
    integer      :: i, j
    call Covariance_Matrix(m,n,X,S)
    call Inverse(m,S,iS)
    !! Forall not used
    do i=1,n
        do j=1,n
            Delta(1,:,j,i) = X(:,i) - X(:,j)
            Temp (1,:,j,i) = Matmul(Delta(1,:,j,i), iS)
            D(i,j) = DSQRT(DOT_PRODUCT(Temp(1,:,j,i), Delta(1,:,j,i)))
        end do
    end do
end subroutine Mahalanobis_Distance

subroutine Cosine(n,u,v,c)
    implicit none
    integer,intent(in)      :: n
    real(kind=8),intent(in) :: u(n), v(n)
    real(kind=8),intent(out):: c
    real(kind=8)    :: l0, l1, l2
    l0 = DOT_PRODUCT(u,v)
    l1 = DOT_PRODUCT(u,u)
    l2 = DOT_PRODUCT(v,v)
    c = l0/(DSQRT(l1)*DSQRT(l2))
end subroutine Cosine

subroutine Hamming_Distance(n,u,v,dist)
    implicit none
    integer,intent(in)      :: n
    real(kind=8),intent(in) :: u(n), v(n)
    real(kind=8),intent(out):: dist
    real(kind=8)            :: diff(n)
    real(kind=8),parameter  :: delta = 1.0D-8
    diff = u - v
    dist = COUNT(DABS(diff)>delta)
end subroutine Hamming_Distance

subroutine Jaccard_Similarity_Coefficient()
    implicit none
    continue
end subroutine Jaccard_Similarity_Coefficient

subroutine Correlation_Coefficient(n,u,v,c)
    implicit none
    integer,intent(in) :: n
    real(kind=8),intent(in) :: u(n)
    real(kind=8),intent(in) :: v(n)
    real(kind=8),intent(out):: c
    real(kind=8) :: su, sv
    real(kind=8) :: denominator, numerator
    su = Standard_Deviation(n,u)
    sv = Standard_Deviation(n,v)
    denominator = su*sv
    numerator = Covariance(n,u,v)
    c = numerator/denominator
end subroutine Correlation_Coefficient

subroutine Tanimoto_Coefficient(n,u,v,t)
    implicit none
    integer,intent(in)      :: n
    real(kind=8),intent(in) :: u(n), v(n)
    real(kind=8),intent(out):: t
    real(kind=8) :: denominator, numerator
    numerator = DOT_PRODUCT(u,v)
    denominator = DOT_PRODUCT(u,u) + DOT_PRODUCT(v,v) - numerator
    t = numerator/denominator
end subroutine Tanimoto_Coefficient

subroutine Lance_Williams_Distance(n,u,v,dist)
    implicit none
    integer,intent(in) :: n
    real(kind=8),intent(in) :: u(n), v(n)
    real(kind=8),intent(out):: dist
    real(kind=8) :: a1(n), a2(n), a3(n)
    a1 = DABS(u-v)
    a2 = u + v
    a3 = a1/a2
    dist = SUM(a3)
end subroutine Lance_Williams_Distance

subroutine Dice_Coefficient()
    !! https://www.biaodianfu.com/dice-coefficient.html
    !! https://baike.baidu.com/item/Dice%E7%B3%BB%E6%95%B0/6597666
    implicit none
    continue
end subroutine Dice_Coefficient

end module Similarity