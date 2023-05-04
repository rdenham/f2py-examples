! Implement the Weber algorithm for a vector median. The functions in 
! this file are intended to be compiled using f2py, and then called from
! the Python module vectormedian.py. These functions are designed to work
! on a stack on multi-band images, which is what defines the shape and structure
! of the input array. Further details on usage can be found in vectormedian.py. 
!
! If no non-null values are found, for a given pixel, then the index is 
! returned as 0 for that pixel. Otherwise, it is an index (starting at 1)
! into the list of input images. It is suitable for use with the selectmedianimage()
! function. 
!
subroutine vectormedian(imagestack, medianNdx, useNull, nullVal, nImgs, nBands, nRows, nCols)
    integer nImgs, nBands, nRows, nCols
    real imagestack(nImgs, nBands, nRows, nCols), nullVal
    integer medianNdx(2, nRows, nCols)
    logical useNull
    ! Directives for f2py
    !f2py intent (in) imagestack, useNull, nullVal
    !f2py intent (out) medianNdx
    !f2py intent (in, hide) nImgs, nBands, nRows, nCols
    real dist(nImgs, nImgs), v1(nBands), v2(nBands), distsum(nImgs)
    integer r, c, i, j, k, ndx(1)
    logical isNull_i, isNull_j, nonNullPts(nImgs)
    
    ! Loop over every row and column
    do r = 1, nRows
        do c = 1, nCols
            ! For a single pixel, across all images in the stack, calculate
            ! the matrix of distances between each image and every other
            ! image, in the nBands-dimensional space. 
            do i = 1, nImgs
                v1 = imagestack(i, :, r, c)
                isNull_i = useNull .and. any(v1.eq.nullVal)
                do j = 1, nImgs
                    v2 = imagestack(j, :, r, c)
                    isNull_j = useNull .and. any(v2.eq.nullVal)
                    if (isNull_i .or. isNull_j) then
                        ! Nulls contribute nothing to the distance sum
                        dist(i, j) = 0
                    elseif (j .gt. i) then
                        d = 0
                        do k = 1, nBands
                            d = d + (v1(k) - v2(k))**2
                        enddo
                        dist(i, j) = sqrt(d)
                    elseif (j .eq. i) then
                        dist(i, j) = 0
                    else
                        ! Save calculation time by using the symmetry of dist array
                        dist(i, j) = dist(j, i)
                    endif
                enddo
                ! Record whether the current point was null
                nonNullPts(i) = .not.isNull_i
            enddo
            
            ! Find the index with the smallest sum of distances to all others
            distsum = sum(dist, 1)
            ndx = minloc(distsum, nonNullPts)
            ! The index is stored as starting from one, suitable for use with selectmedianimage()
            medianNdx(1, r, c) = ndx(1)
            ! A count of how many non-null images we had for this pixel
            medianNdx(2, r, c) = count(nonNullPts)
        enddo
    enddo
end subroutine



! Given the array of indexes returned by the vectormedian() function, use
! these to select out the elements of a given array, creating the actual
! median image. 
! If the medianNdx indicates null (i.e. it is 0), then that pixel of the median
! image will be the null value as given to this function. 
!
subroutine selectmedianimage(imagestack, medianNdx, medianImage, nullVal,                   &
        nImgs, nBands, nRows, nCols)
    integer nImgs, nBands, nRows, nCols
    real imagestack(nImgs, nBands, nRows, nCols), medianImage(nBands, nRows, nCols), nullVal
    real medianNdx(nRows, nCols)
    ! Directives for f2py
    !f2py intent (in) imagestack, medianNdx
    !f2py intent (out) medianImage
    !f2py intent (in, hide) nImgs, nBands, nRows, nCols

    integer r, c, n
    
    do r = 1, nRows
        do c = 1, nCols
            n = medianNdx(r, c)
            if (n .gt. 0) then
                medianImage(:, r, c) = imagestack(n, :, r, c)
            else
                medianImage(:, r, c) = nullVal
            endif
        enddo
    enddo
end subroutine
