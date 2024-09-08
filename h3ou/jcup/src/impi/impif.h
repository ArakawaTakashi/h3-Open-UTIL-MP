! /* -*- fortran -*- */

        integer IMPI_DEFTAG
        integer IMPI_INTEGER
        integer IMPI_INT
        integer IMPI_LONG
        integer IMPI_CHAR
        integer IMPI_FLOAT
        integer IMPI_DOUBLE
        integer IMPI_DOUBLE_PRECISION
        parameter (IMPI_DEFTAG=524288)
        parameter (IMPI_INTEGER=1+524288)
        parameter (IMPI_INT=1+524288)
        parameter (IMPI_LONG=2+524288)
        parameter (IMPI_CHAR=3+524288)
        parameter (IMPI_FLOAT=4+524288)
        parameter (IMPI_DOUBLE=5+524288)
        parameter (IMPI_DOUBLE_PRECISION=5+524288)

        integer IMPI_SUM
        integer IMPI_MAX
        integer IMPI_MIN
        parameter (IMPI_SUM=1+524288)
        parameter (IMPI_MAX=2+524288)
        parameter (IMPI_MIN=3+524288)

        integer WAITIO_STATUS_SIZE
        parameter (WAITIO_STATUS_SIZE=4)
        integer WAITIO_REQUEST_SIZE
        parameter (WAITIO_REQUEST_SIZE=22)
