PROGRAM OpenMP_test2
USE OMP_LIB

INTEGER :: thread_id

!$OMP PARALLEL PRIVATE(thread_id)

    thread_id = OMP_GET_THREAD_NUM()
    PRINT *, "Hello from process: ", thread_id

!$OMP END PARALLEL

END
