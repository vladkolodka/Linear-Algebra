<?php
namespace vladkolodka\linearAlgebra;

class Matrix {
    static function hypo($a, $b) {
        if (abs($a) > abs($b)) {
            $r = $b / $a;
            $r = abs($a) * sqrt(1 + $r * $r);
        } elseif ($b != 0) {
            $r = $a / $b;
            $r = abs($b) * sqrt(1 + $r * $r);
        } else {
            $r = 0.0;
        }
        return $r;
    }


    /**
     * Matrix storage
     *
     * @var array
     * @access private
     */
    public $matrix = array();

    /**
     * @var integer
     */
    public $m;

    /**
     * @var integer
     */
    public $n;

    public function __construct($a = null, $b = null) {
        if (!$a && !$b) throw new MatrixException(0);

        $prototype = implode(',', array_map('gettype', func_get_args()));
        switch ($prototype) {
            // wrap existing array matrix
            case 'array':
                $this->m = count($a);
                $this->n = count($a[0]);

                $this->matrix = $a;
                break;

            // empty matrix
            case 'integer':
                $this->m = $a;
                $this->n = $a;

                $this->matrix = array_fill(0, $a, array_fill(0, $a, 0));
                break;

            // empty matrix
            case 'integer,integer':
                $this->m = $a;
                $this->n = $b;

                $this->matrix = array_fill(0, $a, array_fill(0, $b, 0));
                break;

            // Rectangular matrix
            case 'array,integer':
                $this->m = $b;

                if ($this->m != 0) {
                    $this->n = count($a) / $b;
                } else {
                    $this->n = 0;
                }
                if (($this->m * $this->n) == count($a)) {
                    for ($i = 0; $i < $this->m; ++$i) {
                        for ($j = 0; $j < $this->n; ++$j) {
                            $this->matrix[$i][$j] = $a[$i + $j * $b];
                        }
                    }
                } else {
                    throw new MatrixException(4);
                }
                break;
            default:
                throw new MatrixException(0);
                break;
        }
    }

    public function getArray() {
        return $this->matrix;
    }

    /**
     *    getRowDimension
     *
     * @return integer Row dimension
     */
    public function getRowDimension() {
        return $this->m;
    }

    /**
     *    getColumnDimension
     *
     * @return integer Column dimension
     */
    public function getColumnDimension() {
        return $this->n;
    }

    /**
     * Get value by indexes
     * @param $i
     * @param $j
     * @return integer|float
     */
    public function get($i, $j) {
        return $this->matrix[$i][$j];
    }

    /**
     *    getMatrix
     *
     *    Get a sub-matrix
     * @return Matrix sub-matrix
     * @throws MatrixException
     */
    public function getMatrix() {
        if (func_num_args() > 0) {
            $args = func_get_args();
            $match = implode(",", array_map('gettype', $args));

            switch ($match) {
                //A($i0...; $j0...)
                case 'integer,integer':
                    list($i0, $j0) = $args;
                    if ($i0 >= 0) {
                        $m = $this->m - $i0;
                    } else {
                        throw new MatrixException(2);
                    }
                    if ($j0 >= 0) {
                        $n = $this->n - $j0;
                    } else {
                        throw new MatrixException(2);
                    }
                    $R = new Matrix($m, $n);
                    for ($i = $i0; $i < $this->m; ++$i) {
                        for ($j = $j0; $j < $this->n; ++$j) {
                            $R->set($i, $j, $this->matrix[$i][$j]);
                        }
                    }
                    return $R;
                    break;
                //A($i0...$iF; $j0...$jF)
                case 'integer,integer,integer,integer':
                    list($i0, $iF, $j0, $jF) = $args;
                    if (($iF > $i0) && ($this->m >= $iF) && ($i0 >= 0)) {
                        $m = $iF - $i0;
                    } else {
                        throw new MatrixException(2);
                    }
                    if (($jF > $j0) && ($this->n >= $jF) && ($j0 >= 0)) {
                        $n = $jF - $j0;
                    } else {
                        throw new MatrixException(2);
                    }
                    $R = new Matrix($m + 1, $n + 1);
                    for ($i = $i0; $i <= $iF; ++$i) {
                        for ($j = $j0; $j <= $jF; ++$j) {
                            $R->set($i - $i0, $j - $j0, $this->matrix[$i][$j]);
                        }
                    }
                    return $R;
                    break;
                //$R = array of row indices; $C = array of column indices
                case 'array,array':
                    list($RL, $CL) = $args;
                    if (count($RL) > 0) {
                        $m = count($RL);
                    } else {
                        throw new MatrixException(2);
                    }
                    if (count($CL) > 0) {
                        $n = count($CL);
                    } else {
                        throw new MatrixException(2);
                    }
                    $R = new Matrix($m, $n);
                    for ($i = 0; $i < $m; ++$i) {
                        for ($j = 0; $j < $n; ++$j) {
                            $R->set($i, $j, $this->matrix[$RL[$i]][$CL[$j]]);
                        }
                    }
                    return $R;
                    break;
                //A($i0...$iF); $CL = array of column indices
                case 'integer,integer,array':
                    list($i0, $iF, $CL) = $args;
                    if (($iF > $i0) && ($this->m >= $iF) && ($i0 >= 0)) {
                        $m = $iF - $i0;
                    } else {
                        throw new MatrixException(2);
                    }
                    if (count($CL) > 0) {
                        $n = count($CL);
                    } else {
                        throw new MatrixException(2);
                    }
                    $R = new Matrix($m, $n);
                    for ($i = $i0; $i < $iF; ++$i) {
                        for ($j = 0; $j < $n; ++$j) {
                            $R->set($i - $i0, $j, $this->matrix[$i][$j]);
                        }
                    }
                    return $R;
                    break;
                //$RL = array of row indices
                case 'array,integer,integer':
                    list($RL, $j0, $jF) = $args;
                    if (count($RL) > 0) {
                        $m = count($RL);
                    } else {
                        throw new MatrixException(2);
                    }
                    if (($jF >= $j0) && ($this->n >= $jF) && ($j0 >= 0)) {
                        $n = $jF - $j0;
                    } else {
                        throw new MatrixException(2);
                    }
                    $R = new Matrix($m, $n + 1);
                    for ($i = 0; $i < $m; ++$i) {
                        for ($j = $j0; $j <= $jF; ++$j) {
                            $R->set($i, $j - $j0, $this->matrix[$RL[$i]][$j]);
                        }
                    }
                    return $R;
                    break;
                default:
                    throw new MatrixException(0);
                    break;
            }
        } else {
            throw new MatrixException(0);
        }
    }

    /**
     *    checkMatrixDimensions
     *
     *    Is matrix B the same size?
     * @param Matrix $b Matrix B
     *
     * @throws MatrixException
     * @return boolean
     */
    public function checkMatrixDimensions($b = null) {
        if ($b instanceof Matrix) {
            if (($this->m == $b->getRowDimension()) && ($this->n == $b->getColumnDimension())) {
                return true;
            } else throw new MatrixException(3);

        } else throw new MatrixException(1);
    }

    /**
     *    set
     *
     *    Set the i,j-th element of the matrix.
     * @param int $i Row position
     * @param int $j Column position
     * @param mixed $c Int/float/double value
     * @return mixed Element (int/float/double)
     */
    public function set($i, $j, $c) {
        // Optimized set version just has this
        $this->matrix[$i][$j] = $c;
    }

    /**
     *    identity
     *
     *    Generate an identity matrix.
     * @param int $m Row dimension
     * @param int $n Column dimension
     * @return Matrix Identity matrix
     */
    public function identity($m = null, $n = null) {
        return $this->diagonal($m, $n, 1);
    }

    /**
     *    diagonal
     *
     *    Generate a diagonal matrix
     * @param int $m Row dimension
     * @param int $n Column dimension
     * @param mixed $c Diagonal value
     * @return Matrix Diagonal matrix
     */
    public function diagonal($m = null, $n = null, $c = 1) {
        $R = new Matrix($m, $n);
        for ($i = 0; $i < $m; ++$i) {
            $R->set($i, $i, $c);
        }
        return $R;
    }

    /**
     *    getMatrixByRow
     *
     *    Get a sub-matrix by row index/range
     * @param int $i0 Initial row index
     * @param int $iF Final row index
     *
     * @throws MatrixException
     * @return Matrix Sub-matrix
     */
    public function getMatrixByRow($i0 = null, $iF = null) {
        if (is_int($i0)) {
            if (is_int($iF)) {
                return $this->getMatrix($i0, 0, $iF + 1, $this->n);
            } else {
                return $this->getMatrix($i0, 0, $i0 + 1, $this->n);
            }
        } else {
            throw new MatrixException(1);
        }
    }

    /**
     *    getMatrixByCol
     *
     *    Get a sub-matrix by column index/range
     * @param $j0 int $i0 Initial column index
     * @param $jF int $iF Final column index
     *
     * @throws MatrixException
     * @return Matrix Sub-matrix
     */
    public function getMatrixByCol($j0 = null, $jF = null) {
        if (is_int($j0)) {
            if (is_int($jF)) {
                return $this->getMatrix(0, $j0, $this->m, $jF + 1);
            } else {
                return $this->getMatrix(0, $j0, $this->m, $j0 + 1);
            }
        } else {
            throw new MatrixException(1);
        }
    }

    /**
     *    transpose
     *
     *    Tranpose matrix
     * @return Matrix Transposed matrix
     */
    public function transpose() {
        $R = new Matrix($this->n, $this->m);
        for ($i = 0; $i < $this->m; ++$i) {
            for ($j = 0; $j < $this->n; ++$j) {
                $R->set($j, $i, $this->matrix[$i][$j]);
            }
        }
        return $R;
    }

    /**
     *    trace
     *
     *    Sum of diagonal elements
     * @return float Sum of diagonal elements
     */
    public function trace() {
        $s = 0;
        $n = min($this->m, $this->n);
        for ($i = 0; $i < $n; ++$i) {
            $s += $this->matrix[$i][$i];
        }
        return $s;
    }

    /**
     *    plus
     *
     *    A + B
     * @param mixed $a Matrix/Array
     *
     * @throws MatrixException
     * @return Matrix Sum
     */
    public function plus($a) {
        if (func_num_args() > 0) {
            $args = func_get_args();
            $prototype = implode(",", array_map('gettype', $args));

            switch ($prototype) {
                case 'object':
                    if ($a instanceof Matrix) {
                        $M = $a;
                    } else {
                        throw new MatrixException(1);
                    }
                    break;
                case 'array':
                    $M = new Matrix($a);
                    break;
                default:
                    throw new MatrixException(0);
                    break;
            }
            $this->checkMatrixDimensions($M);
            for ($i = 0; $i < $this->m; ++$i) {
                for ($j = 0; $j < $this->n; ++$j) {
                    $M->set($i, $j, $M->get($i, $j) + $this->matrix[$i][$j]);
                }
            }
            return $M;
        } else throw new MatrixException(0);
    }

    /**
     *    plusEquals
     *
     *    A = A + B
     * @param $a mixed Matrix/Array
     *
     * @throws MatrixException
     * @return Matrix Sum
     */
    public function plusEquals($a) {
        if (func_num_args() > 0) {
            $args = func_get_args();
            $match = implode(",", array_map('gettype', $args));

            switch ($match) {
                case 'object':
                    if ($a instanceof Matrix) {
                        $M = $a;
                    } else {
                        throw new MatrixException(1);
                    }
                    break;
                case 'array':
                    $M = new Matrix($a);
                    break;
                default:
                    throw new MatrixException(1);
                    break;
            }
            $this->checkMatrixDimensions($M);
            for ($i = 0; $i < $this->m; ++$i) {
                for ($j = 0; $j < $this->n; ++$j) {
                    $validValues = true;
                    $value = $M->get($i, $j);
                    if ((is_string($this->matrix[$i][$j])) && (strlen($this->matrix[$i][$j]) > 0) && (!is_numeric($this->matrix[$i][$j]))) {
                        $this->matrix[$i][$j] = trim($this->matrix[$i][$j], '"');
                        $validValues &= $this->matrix[$i][$j];
                    }
                    if ((is_string($value)) && (strlen($value) > 0) && (!is_numeric($value))) {
                        $value = trim($value, '"');
                        $validValues &= $value;
                    }
                    if ($validValues) {
                        $this->matrix[$i][$j] += $value;
                    } else {
                        $this->matrix[$i][$j] = 'NaN';
                    }
                }
            }
            return $this;
        } else {
            throw new MatrixException(0);
        }
    }

    /**
     *    minus
     *
     *    A - B
     * @return Matrix Sum
     * @throws MatrixException
     * @internal param mixed $B Matrix/Array
     */
    public function minus() {
        if (func_num_args() > 0) {
            $args = func_get_args();
            $match = implode(",", array_map('gettype', $args));

            switch ($match) {
                case 'object':
                    if ($args[0] instanceof Matrix) {
                        $M = $args[0];
                    } else {
                        throw new MatrixException(1);
                    }
                    break;
                case 'array':
                    $M = new Matrix($args[0]);
                    break;
                default:
                    throw new MatrixException(0);
                    break;
            }
            $this->checkMatrixDimensions($M);
            for ($i = 0; $i < $this->m; ++$i) {
                for ($j = 0; $j < $this->n; ++$j) {
                    $M->set($i, $j, $M->get($i, $j) - $this->matrix[$i][$j]);
                }
            }
            return $M;
        } else {
            throw new MatrixException(0);
        }
    }

    /**
     *    minusEquals
     *
     *    A = A - B
     * @return Matrix Sum
     * @throws MatrixException
     * @internal param mixed $B Matrix/Array
     */
    public function minusEquals() {
        if (func_num_args() > 0) {
            $args = func_get_args();
            $match = implode(",", array_map('gettype', $args));

            switch ($match) {
                case 'object':
                    if ($args[0] instanceof Matrix) {
                        $M = $args[0];
                    } else {
                        throw new MatrixException(1);
                    }
                    break;
                case 'array':
                    $M = new Matrix($args[0]);
                    break;
                default:
                    throw new MatrixException(0);
                    break;
            }
            $this->checkMatrixDimensions($M);
            for ($i = 0; $i < $this->m; ++$i) {
                for ($j = 0; $j < $this->n; ++$j) {
                    $validValues = true;
                    $value = $M->get($i, $j);
                    if ((is_string($this->matrix[$i][$j])) && (strlen($this->matrix[$i][$j]) > 0) && (!is_numeric($this->matrix[$i][$j]))) {
                        $this->matrix[$i][$j] = trim($this->matrix[$i][$j], '"');
                        $validValues &= $this->matrix[$i][$j];
                    }
                    if ((is_string($value)) && (strlen($value) > 0) && (!is_numeric($value))) {
                        $value = trim($value, '"');
                        $validValues &= $value;
                    }
                    if ($validValues) {
                        $this->matrix[$i][$j] -= $value;
                    } else {
                        $this->matrix[$i][$j] = 'NaN';
                    }
                }
            }
            return $this;
        } else {
            throw new MatrixException(0);
        }
    }

    /**
     *    arrayTimes
     *
     *    Element-by-element multiplication
     *    Cij = Aij * Bij
     * @return Matrix Matrix Cij
     * @throws MatrixException
     * @internal param mixed $B Matrix/Array
     */
    public function arrayTimes() {
        if (func_num_args() > 0) {
            $args = func_get_args();
            $match = implode(",", array_map('gettype', $args));

            switch ($match) {
                case 'object':
                    if ($args[0] instanceof Matrix) {
                        $M = $args[0];
                    } else {
                        throw new MatrixException(1);
                    }
                    break;
                case 'array':
                    $M = new Matrix($args[0]);
                    break;
                default:
                    throw new MatrixException(0);
                    break;
            }
            $this->checkMatrixDimensions($M);
            for ($i = 0; $i < $this->m; ++$i) {
                for ($j = 0; $j < $this->n; ++$j) {
                    $M->set($i, $j, $M->get($i, $j) * $this->matrix[$i][$j]);
                }
            }
            return $M;
        } else {
            throw new MatrixException(0);
        }
    }

    /**
     *    arrayTimesEquals
     *
     *    Element-by-element multiplication
     *    Aij = Aij * Bij
     * @return Matrix Matrix Aij
     * @throws MatrixException
     * @internal param mixed $B Matrix/Array
     */
    public function arrayTimesEquals() {
        if (func_num_args() > 0) {
            $args = func_get_args();
            $match = implode(",", array_map('gettype', $args));

            switch ($match) {
                case 'object':
                    if ($args[0] instanceof Matrix) {
                        $M = $args[0];
                    } else {
                        throw new MatrixException(1);
                    }
                    break;
                case 'array':
                    $M = new Matrix($args[0]);
                    break;
                default:
                    throw new MatrixException(0);
                    break;
            }
            $this->checkMatrixDimensions($M);
            for ($i = 0; $i < $this->m; ++$i) {
                for ($j = 0; $j < $this->n; ++$j) {
                    $validValues = true;
                    $value = $M->get($i, $j);
                    if ((is_string($this->matrix[$i][$j])) && (strlen($this->matrix[$i][$j]) > 0) && (!is_numeric($this->matrix[$i][$j]))) {
                        $this->matrix[$i][$j] = trim($this->matrix[$i][$j], '"');
                        $validValues &= $this->matrix[$i][$j];
                    }
                    if ((is_string($value)) && (strlen($value) > 0) && (!is_numeric($value))) {
                        $value = trim($value, '"');
                        $validValues &= $value;
                    }
                    if ($validValues) {
                        $this->matrix[$i][$j] *= $value;
                    } else {
                        $this->matrix[$i][$j] = 'NaN';
                    }
                }
            }
            return $this;
        } else {
            throw new MatrixException(0);
        }
    }

    /**
     *    arrayRightDivide
     *
     *    Element-by-element right division
     *    A / B
     * @return Matrix Division result
     * @throws MatrixException
     * @internal param Matrix $B Matrix B
     */
    public function arrayRightDivide() {
        if (func_num_args() > 0) {
            $args = func_get_args();
            $match = implode(",", array_map('gettype', $args));

            switch ($match) {
                case 'object':
                    if ($args[0] instanceof Matrix) {
                        $M = $args[0];
                    } else {
                        throw new MatrixException(1);
                    }
                    break;
                case 'array':
                    $M = new Matrix($args[0]);
                    break;
                default:
                    throw new MatrixException(0);
                    break;
            }
            $this->checkMatrixDimensions($M);
            for ($i = 0; $i < $this->m; ++$i) {
                for ($j = 0; $j < $this->n; ++$j) {
                    $validValues = true;
                    $value = $M->get($i, $j);
                    if ((is_string($this->matrix[$i][$j])) && (strlen($this->matrix[$i][$j]) > 0) && (!is_numeric($this->matrix[$i][$j]))) {
                        $this->matrix[$i][$j] = trim($this->matrix[$i][$j], '"');
                        $validValues &= $this->matrix[$i][$j];
                    }
                    if ((is_string($value)) && (strlen($value) > 0) && (!is_numeric($value))) {
                        $value = trim($value, '"');
                        $validValues &= $value;
                    }
                    if ($validValues) {
                        if ($value == 0) {
                            //    Trap for Divide by Zero error
                            $M->set($i, $j, '#DIV/0!');
                        } else {
                            $M->set($i, $j, $this->matrix[$i][$j] / $value);
                        }
                    } else {
                        $M->set($i, $j, 'MaM');
                    }
                }
            }
            return $M;
        } else {
            throw new MatrixException(0);
        }
    }


    /**
     *    arrayRightDivideEquals
     *
     *    Element-by-element right division
     *    Aij = Aij / Bij
     * @return Matrix Matrix Aij
     * @throws MatrixException
     * @internal param mixed $B Matrix/Array
     */
    public function arrayRightDivideEquals() {
        if (func_num_args() > 0) {
            $args = func_get_args();
            $match = implode(",", array_map('gettype', $args));

            switch ($match) {
                case 'object':
                    if ($args[0] instanceof Matrix) {
                        $M = $args[0];
                    } else {
                        throw new MatrixException(1);
                    }
                    break;
                case 'array':
                    $M = new Matrix($args[0]);
                    break;
                default:
                    throw new MatrixException(0);
                    break;
            }
            $this->checkMatrixDimensions($M);
            for ($i = 0; $i < $this->m; ++$i) {
                for ($j = 0; $j < $this->n; ++$j) {
                    $this->matrix[$i][$j] = $this->matrix[$i][$j] / $M->get($i, $j);
                }
            }
            return $M;
        } else {
            throw new MatrixException(0);
        }
    }


    /**
     *    arrayLeftDivide
     *
     *    Element-by-element Left division
     *    A / B
     * @return Matrix Division result
     * @throws MatrixException
     * @internal param Matrix $B Matrix B
     */
    public function arrayLeftDivide() {
        if (func_num_args() > 0) {
            $args = func_get_args();
            $match = implode(",", array_map('gettype', $args));

            switch ($match) {
                case 'object':
                    if ($args[0] instanceof Matrix) {
                        $M = $args[0];
                    } else {
                        throw new MatrixException(1);
                    }
                    break;
                case 'array':
                    $M = new Matrix($args[0]);
                    break;
                default:
                    throw new MatrixException(0);
                    break;
            }
            $this->checkMatrixDimensions($M);
            for ($i = 0; $i < $this->m; ++$i) {
                for ($j = 0; $j < $this->n; ++$j) {
                    $M->set($i, $j, $M->get($i, $j) / $this->matrix[$i][$j]);
                }
            }
            return $M;
        } else {
            throw new MatrixException(0);
        }
    }


    /**
     *    arrayLeftDivideEquals
     *
     *    Element-by-element Left division
     *    Aij = Aij / Bij
     * @return Matrix Matrix Aij
     * @throws MatrixException
     * @internal param mixed $B Matrix/Array
     */
    public function arrayLeftDivideEquals() {
        if (func_num_args() > 0) {
            $args = func_get_args();
            $match = implode(",", array_map('gettype', $args));

            switch ($match) {
                case 'object':
                    if ($args[0] instanceof Matrix) {
                        $M = $args[0];
                    } else {
                        throw new MatrixException(1);
                    }
                    break;
                case 'array':
                    $M = new Matrix($args[0]);
                    break;
                default:
                    throw new MatrixException(0);
                    break;
            }
            $this->checkMatrixDimensions($M);
            for ($i = 0; $i < $this->m; ++$i) {
                for ($j = 0; $j < $this->n; ++$j) {
                    $this->matrix[$i][$j] = $M->get($i, $j) / $this->matrix[$i][$j];
                }
            }
            return $M;
        } else {
            throw new MatrixException(0);
        }
    }


    /**
     *    times
     *
     *    Matrix multiplication
     * @return Matrix Product
     * @throws MatrixException
     * @internal param mixed $n Matrix/Array/Scalar
     */
    public function times() {
        if (func_num_args() > 0) {
            $args = func_get_args();
            $prototype = implode(",", array_map('gettype', $args));

            $B = null;
            switch ($prototype) {
                case 'object':
                    if ($args[0] instanceof Matrix) {
                        $B = $args[0];
                    } else {
                        throw new MatrixException(1);
                    }
                    if ($this->n == $B->m) {
                        $C = new Matrix($this->m, $B->n);
                        $bColJ = [];
                        for ($j = 0; $j < $B->n; ++$j) {
                            for ($k = 0; $k < $this->n; ++$k) {
                                $bColJ[$k] = $B->matrix[$k][$j];
                            }
                            for ($i = 0; $i < $this->m; ++$i) {
                                $aRowI = $this->matrix[$i];
                                $s = 0;
                                for ($k = 0; $k < $this->n; ++$k) {
                                    $s += $aRowI[$k] * $bColJ[$k];
                                }
                                $C->matrix[$i][$j] = $s;
                            }
                        }
                        return $C;
                    } else {
                        throw new MatrixException(3);
                    }
                    break;
                case 'array':
                    $B = new Matrix($args[0]);
                    if ($this->n == $B->m) {
                        $C = new Matrix($this->m, $B->n);
                        for ($i = 0; $i < $C->m; ++$i) {
                            for ($j = 0; $j < $C->n; ++$j) {
                                $s = "0";
                                for ($k = 0; $k < $C->n; ++$k) {
                                    $s += $this->matrix[$i][$k] * $B->matrix[$k][$j];
                                }
                                $C->matrix[$i][$j] = $s;
                            }
                        }
                        return $C;
                    } else {
                        throw new MatrixException(3);
                    }
                    break;
                case 'integer':
                    $C = new Matrix($this->matrix);
                    for ($i = 0; $i < $C->m; ++$i) {
                        for ($j = 0; $j < $C->n; ++$j) {
                            $C->matrix[$i][$j] *= $args[0];
                        }
                    }
                    return $C;
                    break;
                case 'double':
                    $C = new Matrix($this->m, $this->n);
                    for ($i = 0; $i < $C->m; ++$i) {
                        for ($j = 0; $j < $C->n; ++$j) {
                            $C->matrix[$i][$j] = $args[0] * $this->matrix[$i][$j];
                        }
                    }
                    return $C;
                    break;
                case 'float':
                    $C = new Matrix($this->matrix);
                    for ($i = 0; $i < $C->m; ++$i) {
                        for ($j = 0; $j < $C->n; ++$j) {
                            $C->matrix[$i][$j] *= $args[0];
                        }
                    }
                    return $C;
                    break;
                default:
                    throw new MatrixException(0);
                    break;
            }
        } else {
            throw new MatrixException(0);
        }
    }

    /**
     *    power
     *
     *    A = A ^ B
     * @return Matrix Sum
     * @throws MatrixException
     * @internal param mixed $B Matrix/Array
     */
    public function power() {
        if (func_num_args() > 0) {
            $args = func_get_args();
            $match = implode(",", array_map('gettype', $args));

            switch ($match) {
                case 'object':
                    if ($args[0] instanceof Matrix) {
                        $M = $args[0];
                    } else {
                        throw new MatrixException(1);
                    }
                    break;
                case 'array':
                    $M = new Matrix($args[0]);
                    break;
                default:
                    throw new MatrixException(0);
                    break;
            }
            $this->checkMatrixDimensions($M);
            for ($i = 0; $i < $this->m; ++$i) {
                for ($j = 0; $j < $this->n; ++$j) {
                    $validValues = true;
                    $value = $M->get($i, $j);
                    if ((is_string($this->matrix[$i][$j])) && (strlen($this->matrix[$i][$j]) > 0) && (!is_numeric($this->matrix[$i][$j]))) {
                        $this->matrix[$i][$j] = trim($this->matrix[$i][$j], '"');
                        $validValues &= $this->matrix[$i][$j];
                    }
                    if ((is_string($value)) && (strlen($value) > 0) && (!is_numeric($value))) {
                        $value = trim($value, '"');
                        $validValues &= $value;
                    }
                    if ($validValues) {
                        $this->matrix[$i][$j] = pow($this->matrix[$i][$j], $value);
                    } else {
                        $this->matrix[$i][$j] = 'NaN';
                    }
                }
            }
            return $this;
        } else {
            throw new MatrixException(0);
        }
    }

    /**
     *    concat
     *
     *    A = A & B
     * @param $a
     * @return Matrix Sum
     * @throws MatrixException
     * @internal param mixed $B Matrix/Array
     */
    public function concat($a) {
        if (func_num_args() > 0) {
            $prototype = implode(",", array_map('gettype', func_get_args()));

            switch ($prototype) {
                case 'object':
                    if ($a instanceof Matrix) {
                        $M = $a;
                    } else {
                        throw new MatrixException(1);
                    }
                break;
                case 'array':
                    $M = new Matrix($a);
                    break;
                default:
                    throw new MatrixException(0);
                    break;
            }
            $this->checkMatrixDimensions($M);
            for ($i = 0; $i < $this->m; ++$i) {
                for ($j = 0; $j < $this->n; ++$j) {
                    $this->matrix[$i][$j] = trim($this->matrix[$i][$j], '"') . trim($M->get($i, $j), '"');
                }
            }
            return $this;
        } else {
            throw new MatrixException(0);
        }
    }

    /**
     *    Solve A*X = B.
     *
     * @param Matrix $B Right hand side
     * @return Matrix ... Solution if A is square, least squares solution otherwise
     */
    public function solve($B) {
        if ($this->m == $this->n) {
            $LU = new UDecomposition($this);
            return $LU->solve($B);
        } else {
            $QR = new QRDecomposition($this);
            return $QR->solve($B);
        }
    }

    /**
     *    Matrix inverse or pseudo inverse.
     *
     * @return Matrix ... Inverse(A) if A is square, pseudo inverse otherwise.
     */
    public function inverse() {
        return $this->solve($this->identity($this->m, $this->m));
    }

    /**
     *    det
     *
     *    Calculate determinant
     * @return float Determinant
     */
    public function det() {
        $L = new UDecomposition($this);
        return $L->det();
    }
}