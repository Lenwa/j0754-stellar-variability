"""Prevent empty/error responses and metadata wrappers becoming detections."""
import unittest
import query_multiwavelength_archives as q

class ResponseTests(unittest.TestCase):
    def test_mast_and_photon_rows_have_same_meaning(self):
        self.assertEqual(q.json_row_count({'data':[{},{}]}),2)
        self.assertEqual(q.json_row_count({'data':{'Tables':[{'Rows':[[1],[2],[3]]}]}}),3)

    def test_sql_query_receipt_is_not_an_observed_star(self):
        self.assertEqual(q.json_row_count([{'TableName':'Table1','Rows':[]},
            {'TableName':'SqlQuery','Rows':[{'query':'SELECT...'}]}]),0)

    def test_non_tabular_limit_is_not_zero_detections(self):
        self.assertIsNone(q.json_row_count({'Exposure':243.,'Flag_pos':1}))

if __name__=='__main__':unittest.main()
