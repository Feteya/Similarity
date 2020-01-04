!! Test data from
!! https://www.cnblogs.com/liujinhong/p/6001997.html
!! https://www.jianshu.com/p/d2c11c898f72
program Test
    use Similarity
    implicit none

    !! Test Euclidean_Distance, Manhattan_Distance, Chebyshev_Distance, Minkowski_Distance
    real(kind=8)    :: a(2),b(2),c(2)
    !! Test Standardized_Euclidean_Distance
    real(kind=8)    :: s(2)
    real(kind=8)    :: dist

    !! Test Covariance_Matrix
    real(kind=8)    :: M(3,10)
    real(kind=8)    :: COV(3,3)

    !! Test Inverse
    real(kind=8)    :: T(3,3), iT(3,3)

    !! Test 
    real(kind=8)    :: X(2,4)
    real(kind=8)    :: D(4,4)

    !! Test Cosine
    real(kind=8) :: p1(2),p2(2),p3(2)
    real(kind=8) :: cos_p(3)

    !! Test Hamming_Distance
    real(kind=8) :: HD1(2),HD2(2),HD3(2)
    real(kind=8) :: HD(3)

    !! Correlation_Coefficient
    real(kind=8) :: CCX(4), CCY(4),CC

    integer :: i

    write(*,"(A)") "Test begin!"
    write(*,*)

    a = [0.0d0, 0.0d0]
    b = [1.0d0, 0.0d0]
    c = [0.0d0, 2.0d0]
    write(*,"(A,2(1X,F6.3),A)") "a = [", a, " ]"
    write(*,"(A,2(1X,F6.3),A)") "b = [", b, " ]"
    write(*,"(A,2(1X,F6.3),A)") "c = [", c, " ]"
    write(*,*)


    write(*,"(A)") "Euclidean_Distance"
    call Euclidean_Distance(2,a,b,dist)
    write(*,"(A,F0.3)") "dist(a,b) = ",dist
    call Euclidean_Distance(2,a,c,dist)
    write(*,"(A,F0.3)") "dist(a,c) = ",dist
    call Euclidean_Distance(2,b,c,dist)
    write(*,"(A,F0.3)") "dist(b,c) = ",dist
    write(*,*)

    write(*,"(A)") "Manhattan_Distance"
    call Manhattan_Distance(2,a,b,dist)
    write(*,"(A,F0.3)") "dist(a,b) = ",dist
    call Manhattan_Distance(2,a,c,dist)
    write(*,"(A,F0.3)") "dist(a,c) = ",dist
    call Manhattan_Distance(2,b,c,dist)
    write(*,"(A,F0.3)") "dist(b,c) = ",dist
    write(*,*)

    write(*,"(A)") "Chebyshev_Distance"
    call Chebyshev_Distance(2,a,b,dist)
    write(*,"(A,F0.3)") "dist(a,b) = ",dist
    call Chebyshev_Distance(2,a,c,dist)
    write(*,"(A,F0.3)") "dist(a,c) = ",dist
    call Chebyshev_Distance(2,b,c,dist)
    write(*,"(A,F0.3)") "dist(b,c) = ",dist
    write(*,*)

    write(*,"(A)") "Minkowski_Distance, p=1, Manhattan_Distance"
    call Minkowski_Distance(2,a,b,dist,1)
    write(*,"(A,F0.3)") "dist(a,b) = ",dist
    call Minkowski_Distance(2,a,c,dist,1)
    write(*,"(A,F0.3)") "dist(a,c) = ",dist
    call Minkowski_Distance(2,b,c,dist,1)
    write(*,"(A,F0.3)") "dist(b,c) = ",dist
    write(*,*)

    write(*,"(A)") "Minkowski_Distance, p=2, Euclidean_Distance"
    call Minkowski_Distance(2,a,b,dist,2)
    write(*,"(A,F0.3)") "dist(a,b) = ",dist
    call Minkowski_Distance(2,a,c,dist,2)
    write(*,"(A,F0.3)") "dist(a,c) = ",dist
    call Minkowski_Distance(2,b,c,dist,2)
    write(*,"(A,F0.3)") "dist(b,c) = ",dist
    write(*,*)

    write(*,"(A)") "Minkowski_Distance, p=3, Chebyshev_Distance"
    call Minkowski_Distance(2,a,b,dist,3)
    write(*,"(A,F0.3)") "dist(a,b) = ",dist
    call Minkowski_Distance(2,a,c,dist,3)
    write(*,"(A,F0.3)") "dist(a,c) = ",dist
    call Minkowski_Distance(2,b,c,dist,3)
    write(*,"(A,F0.3)") "dist(b,c) = ",dist
    write(*,*)

    s = [0.5d0, 1.0d0]
    write(*,"(A)") "Standardized_Euclidean_Distance, s=[0.5,1.0]"
    call Standardized_Euclidean_Distance(2,a,b,dist,s)
    write(*,"(A,F0.3)") "dist(a,b) = ",dist
    call Standardized_Euclidean_Distance(2,a,c,dist,s)
    write(*,"(A,F0.3)") "dist(a,c) = ",dist
    call Standardized_Euclidean_Distance(2,b,c,dist,s)
    write(*,"(A,F0.3)") "dist(b,c) = ",dist
    write(*,*)

    M(:,1)  = [49.0d0,  7.0d0, 29.0d0]
    M(:,2)  = [ 8.0d0, 19.0d0, 16.0d0]
    M(:,3)  = [12.0d0,  8.0d0, 14.0d0]
    M(:,4)  = [19.0d0, 37.0d0, 22.0d0]
    M(:,5)  = [ 3.0d0, 43.0d0, 21.0d0]
    M(:,6)  = [34.0d0, 17.0d0, 17.0d0]
    M(:,7)  = [20.0d0, 34.0d0, 27.0d0]
    M(:,8)  = [49.0d0, 14.0d0, 37.0d0]
    M(:,9)  = [20.0d0, 26.0d0, 21.0d0]
    M(:,10) = [31.0d0, 41.0d0, 21.0d0]
    write(*,"(A)") "Covariance_Matrix"
    call Covariance_Matrix(3,10,M,COV)
    do i=1,3
        write(*,"(3(F0.4, 4X))") COV(:,i)
    end do
    write(*,*)

    T(:,1) = [1.0d0, 1.0d0, 1.0d0]
    T(:,2) = [1.0d0, 2.0d0, 5.0d0]
    T(:,3) = [1.0d0, 3.0d0, 1.0d0]
    write(*,"(A)") "Inverse"
    call Inverse(3,T,iT)
    do i=1,3
        write(*,"(3(F0.4, 4X))") iT(i,:)
    end do
    write(*,*)


    X(:,1) = [1.0d0, 2.0d0]
    X(:,2) = [1.0d0, 3.0d0]
    X(:,3) = [2.0d0, 2.0d0]
    X(:,4) = [3.0d0, 1.0d0]
    write(*,"(A)") "Mahalanobis_Distance"
    call Mahalanobis_Distance(2,4,X,D)
    write(*,"(A,F0.3)") "dist(1,2) = ",D(1,2)
    write(*,"(A,F0.3)") "dist(1,3) = ",D(1,3)
    write(*,"(A,F0.3)") "dist(1,4) = ",D(1,4)
    write(*,"(A,F0.3)") "dist(2,3) = ",D(2,3)
    write(*,"(A,F0.3)") "dist(2,4) = ",D(2,4)
    write(*,"(A,F0.3)") "dist(3,4) = ",D(3,4)
    write(*,*)

    p1 = [1.0d0,  0.0d0]
    p2 = [1.0d0,  1.732d0]
    p3 = [-1.0d0, 0.0d0]
    write(*,"(A)") "Cosine"
    call Cosine(2,p1,p2,cos_p(1))
    call Cosine(2,p1,p3,cos_p(2))
    call Cosine(2,p2,p3,cos_p(3))
    write(*,"(A,F0.3)") "cos(1,2) = ",cos_p(1)
    write(*,"(A,F0.3)") "cos(1,3) = ",cos_p(2)
    write(*,"(A,F0.3)") "cos(2,3) = ",cos_p(3)
    write(*,*)

    HD1 = [0.0d0, 0.0d0]
    HD2 = [1.0d0, 0.0d0]
    HD3 = [0.0d0, 2.0d0]
    write(*,"(A)") "Hamming_Distance"
    call Hamming_Distance(2,HD1,HD2,HD(1))
    call Hamming_Distance(2,HD1,HD3,HD(2))
    call Hamming_Distance(2,HD2,HD3,HD(3))
    write(*,"(A,F0.3)") "dist(1,2) = ",HD(1)
    write(*,"(A,F0.3)") "dist(1,3) = ",HD(2)
    write(*,"(A,F0.3)") "dist(2,3) = ",HD(3)
    write(*,*)

    CCX = [1.0D0, 2.0D0, 3.0D0, 4.0D0]
    CCY = [3.0D0, 8.0D0, 7.0D0, 6.0D0]
    write(*,"(A)") "Correlation_Coefficient"
    call Correlation_Coefficient(4,CCX,CCY,CC)
    write(*,"(A,F0.3)") "Correlation_Coefficient = ",CC
    write(*,*)


end program Test