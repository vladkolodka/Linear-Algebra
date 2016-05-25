<?php
namespace vladkolodka\linearAlgebra;

class CholeskyDecomposition {
    /**
     *    Decomposition storage
     * @var array
     * @access private
     */
    private $L = array();

    /**
     *    Matrix row and column dimension
     * @var int
     * @access private
     */
    private $m;

    /**
     *    Symmetric positive definite flag
     * @var boolean
     * @access private
     */
    private $isSpd = true;

    /**
     *    CholeskyDecomposition
     *
     *    Class constructor - decomposes symmetric positive definite matrix
     * @param mixed Matrix square symmetric positive definite matrix
     * @throws MatrixException
     */
    public function __construct($A = null) {
        if ($A instanceof Matrix) {
            $this->L = $A->getArray();
            $this->m = $A->getRowDimension();

            for ($i = 0; $i < $this->m; ++$i) {
                for ($j = $i; $j < $this->m; ++$j) {
                    for ($sum = $this->L[$i][$j], $k = $i - 1; $k >= 0; --$k) {
                        $sum -= $this->L[$i][$k] * $this->L[$j][$k];
                    }
                    if ($i == $j) {
                        if ($sum >= 0) {
                            $this->L[$i][$i] = sqrt($sum);
                        } else {
                            $this->isSpd = false;
                        }
                    } else {
                        if ($this->L[$i][$i] != 0) {
                            $this->L[$j][$i] = $sum / $this->L[$i][$i];
                        }
                    }
                }

                for ($k = $i + 1; $k < $this->m; ++$k) {
                    $this->L[$i][$k] = 0.0;
                }
            }
        } else {
            throw new MatrixException(1);
        }
    }

    /**
     *    Is the matrix symmetric and positive definite?
     *
     * @return boolean
     */
    public function isSPD() {
        return $this->isSpd;
    }    //    function isSPD()

    /**
     *    getL
     *
     *    Return triangular factor.
     * @return Matrix Lower triangular matrix
     */
    public function getL() {
        return new Matrix($this->L);
    }    //    function getL()

    /**
     *    Solve A*X = B
     *
     * @param $B Matrix Row-equal matrix
     * @return Matrix L * L' * X = B
     * @throws MatrixException
     */
    public function solve($B = null) {
        if ($B instanceof Matrix) {
            if ($B->getRowDimension() == $this->m) {
                if ($this->isSpd) {
                    $X = $B->getArray();
                    $nx = $B->getColumnDimension();

                    for ($k = 0; $k < $this->m; ++$k) {
                        for ($i = $k + 1; $i < $this->m; ++$i) {
                            for ($j = 0; $j < $nx; ++$j) {
                                $X[$i][$j] -= $X[$k][$j] * $this->L[$i][$k];
                            }
                        }
                        for ($j = 0; $j < $nx; ++$j) {
                            $X[$k][$j] /= $this->L[$k][$k];
                        }
                    }

                    for ($k = $this->m - 1; $k >= 0; --$k) {
                        for ($j = 0; $j < $nx; ++$j) {
                            $X[$k][$j] /= $this->L[$k][$k];
                        }
                        for ($i = 0; $i < $k; ++$i) {
                            for ($j = 0; $j < $nx; ++$j) {
                                $X[$i][$j] -= $X[$k][$j] * $this->L[$k][$i];
                            }
                        }
                    }

                    return new Matrix($X, $this->m, $nx);
                } else throw new MatrixException(8);
            } else throw new MatrixException(3);
        } else throw new MatrixException(1);
    }    //    function solve()
}