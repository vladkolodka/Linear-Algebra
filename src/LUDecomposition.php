<?php
namespace vladkolodka\linearAlgebra;

use vladkolodka\linearAlgebra\Exceptions\MatrixException;

class UDecomposition {
    /**
     *    Decomposition storage
     * @var array
     */
    private $LU = array();

    /**
     *    Row dimension.
     * @var int
     */
    private $m;

    /**
     *    Column dimension.
     * @var int
     */
    private $n;

    /**
     *    Pivot sign.
     * @var int
     */
    private $PiVSign;

    /**
     *    Internal storage of pivot vector.
     * @var array
     */
    private $piv = array();

    /**
     *    LU Decomposition constructor.
     *
     * @throws MatrixException
     * @param $A Matrix
     */
    public function __construct($A) {
        if ($A instanceof Matrix) {
            // Use a "left-looking", dot-product, Crout/Doolittle algorithm.
            $this->LU = $A->getArray();
            $this->m = $A->getRowDimension();
            $this->n = $A->getColumnDimension();
            for ($i = 0; $i < $this->m; ++$i) {
                $this->piv[$i] = $i;
            }
            $this->PiVSign = 1;
            $LU_colJ = array();

            // Outer loop.
            for ($j = 0; $j < $this->n; ++$j) {
                // Make a copy of the j-th column to localize references.
                for ($i = 0; $i < $this->m; ++$i) {
                    $LU_colJ[$i] = &$this->LU[$i][$j];
                }
                // Apply previous transformations.
                for ($i = 0; $i < $this->m; ++$i) {
                    $LU_rowI = $this->LU[$i];
                    // Most of the time is spent in the following dot product.
                    $kMax = min($i, $j);
                    $s = 0.0;
                    for ($k = 0; $k < $kMax; ++$k) {
                        $s += $LU_rowI[$k] * $LU_colJ[$k];
                    }
                    $LU_rowI[$j] = $LU_colJ[$i] -= $s;
                }
                // Find pivot and exchange if necessary.
                $p = $j;
                for ($i = $j + 1; $i < $this->m; ++$i) {
                    if (abs($LU_colJ[$i]) > abs($LU_colJ[$p])) {
                        $p = $i;
                    }
                }
                if ($p != $j) {
                    for ($k = 0; $k < $this->n; ++$k) {
                        $t = $this->LU[$p][$k];
                        $this->LU[$p][$k] = $this->LU[$j][$k];
                        $this->LU[$j][$k] = $t;
                    }
                    $k = $this->piv[$p];
                    $this->piv[$p] = $this->piv[$j];
                    $this->piv[$j] = $k;
                    $this->PiVSign = $this->PiVSign * -1;
                }
                // Compute multipliers.
                if (($j < $this->m) && ($this->LU[$j][$j] != 0.0)) {
                    for ($i = $j + 1; $i < $this->m; ++$i) {
                        $this->LU[$i][$j] /= $this->LU[$j][$j];
                    }
                }
            }
        } else {
            throw new MatrixException(1);
        }
    }

    /**
     *    Get lower triangular factor.
     *
     * @return array Lower triangular factor
     */
    public function getL() {
        $L = [];
        for ($i = 0; $i < $this->m; ++$i) {
            for ($j = 0; $j < $this->n; ++$j) {
                if ($i > $j) {
                    $L[$i][$j] = $this->LU[$i][$j];
                } elseif ($i == $j) {
                    $L[$i][$j] = 1.0;
                } else {
                    $L[$i][$j] = 0.0;
                }
            }
        }
        return new Matrix($L);
    }

    /**
     *    Get upper triangular factor.
     *
     * @return array Upper triangular factor
     */
    public function getU() {
        $U = [];
        for ($i = 0; $i < $this->n; ++$i) {
            for ($j = 0; $j < $this->n; ++$j) {
                if ($i <= $j) {
                    $U[$i][$j] = $this->LU[$i][$j];
                } else {
                    $U[$i][$j] = 0.0;
                }
            }
        }
        return new Matrix($U);
    }

    /**
     *    Return pivot permutation vector.
     *
     * @return array Pivot vector
     */
    public function getPivot() {
        return $this->piv;
    }

    /**
     *    Alias for getPivot
     *
     * @see getPivot
     */
    public function getDoublePivot() {
        return $this->getPivot();
    }    //    function getDoublePivot()

    /**
     *    Is the matrix non-singular?
     *
     * @return true if U, and hence A, is non-singular.
     */
    public function isNonSingular() {
        for ($j = 0; $j < $this->n; ++$j) {
            if ($this->LU[$j][$j] == 0) {
                return false;
            }
        }
        return true;
    }

    /**
     *    Count determinants
     *
     * @throws MatrixException
     * @return array d matrix determinant
     */
    public function det() {
        if ($this->m == $this->n) {
            $d = $this->PiVSign;
            for ($j = 0; $j < $this->n; ++$j) {
                $d *= $this->LU[$j][$j];
            }
            return $d;
        } else {
            throw new MatrixException(3);
        }
    }    //    function det()

    /**
     *    Solve A*X = B
     *
     * @param  $B  Matrix with as many rows as A and any number of columns.
     * @return integer X so that L*U*X = B(piv,:)
     * @throws MatrixException  IllegalArgumentException Matrix row dimensions must agree.
     * @throws MatrixException  RuntimeException  Matrix is singular.
     */
    public function solve($B) {
        if ($B->getRowDimension() == $this->m) {
            if ($this->isNonSingular()) {
                // Copy right hand side with pivoting
                $nx = $B->getColumnDimension();
                $X = $B->getMatrix($this->piv, 0, $nx - 1)->getArray();
                // Solve L*Y = B(piv,:)
                for ($k = 0; $k < $this->n; ++$k) {
                    for ($i = $k + 1; $i < $this->n; ++$i) {
                        for ($j = 0; $j < $nx; ++$j) {
                            $X[$i][$j] -= $X[$k][$j] * $this->LU[$i][$k];
                        }
                    }
                }
                // Solve U*X = Y;
                for ($k = $this->n - 1; $k >= 0; --$k) {
                    for ($j = 0; $j < $nx; ++$j) {
                        $X[$k][$j] /= $this->LU[$k][$k];
                    }
                    for ($i = 0; $i < $k; ++$i) {
                        for ($j = 0; $j < $nx; ++$j) {
                            $X[$i][$j] -= $X[$k][$j] * $this->LU[$i][$k];
                        }
                    }
                }
                return new Matrix($X);
            } else throw new MatrixException(6);
        } else throw new MatrixException(7);
    }
}
