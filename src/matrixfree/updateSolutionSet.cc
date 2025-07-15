#ifndef UPDATE_SOLUTIONSET_H
#define UPDATE_SOLUTIONSET_H

template<int dim>
void MatrixFreePDE<dim> :: updateSolutionSet(std::vector< vectorType * >& dst, std::vector< vectorType * >& src, std::string flag)
{
  for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
    for (unsigned int dof=0; dof<src[fieldIndex]->local_size(); ++dof){
      dst[fieldIndex]->local_element(dof) = src[fieldIndex]->local_element(dof);
      //std::cout << flag << " field " << fieldIndex << "  , solution : " << src[fieldIndex]->local_element(dof) << std::endl;
      //std::cout << flag << " field " << fieldIndex << "  , old_solution : " << dst[fieldIndex]->local_element(dof) << std::endl;
    }
  }
}

#endif
