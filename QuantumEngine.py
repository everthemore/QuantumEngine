import numpy as np
import scipy as sp
import scipy.sparse

class InvalidState(Exception):
    pass

def convertToMixedBase(i, bases):
    representation = []

    if( i >= np.product(bases) ):
        raise InvalidState("This state (%d) can not be represented in this basis"%i)

    for b in bases:
        representation.append( i % b )
        i //= b

    reverse = representation[::-1]
    reversestring = [str(a) for a in reverse]
    return "".join(reversestring)

class HilbertSpace:

    def __init__( self, localHilbertSizes ):
        # Initializes the Hilbert space, defined by the list localHilbertSizes
        # which is a list of integers describing the local Hilbert space dimensions
        # for each qudit

        # Store the local hilbert space size info
        self.localHilbertSizes = localHilbertSizes
        self.N = np.product(localHilbertSizes)

        # Build fast indexing dictionariesself.state_to_index = {}
        self.state_to_index = {}
        self.index_to_state = {}
        for i in range(self.size):
            # Convert this integer to the respective bases
            mixedBaseRepresentation = convertToMixedBase(i, self.localHilbertSizes)
            self.state_to_index[mixedBaseRepresentation] = i;
            self.index_to_state[i] = mixedBaseRepresentation;

        # Initialize the current state to the 'all 0's' state
        self.state = sp.sparse.lil_matrix( (self.N,1), dtype=np.complex_ )
        self.state[ self.state_to_index["0"*len(localHilbertSizes)] ] = 1 + 0j

    def setProductState( self, state ):
        # Sets the Hilbert space to a product state, where state
        # is a list of integers.

        # TODO: Check that state is a list (or make override function for strings)

        # Assert that this is a valid state for the Hilbert space
        for i,s in enumerate(state):
          if s > self.localHilbertSizes[i]:
            raise InvalidState("The state %s is not part of this Hilbert space"%("".join(state)))

#        self.amplitudes = {}
#        self.amplitudes["".join(state)] = 1.0 + 0.0j

        self.state = sp.sparse.lil_matrix( (self.N,1), dtype=np.complex_ )
        self.state[ self.state_to_index["".join([str(s) for s in state])] ] = 1 + 0j

    def constructImaginarySwap(self, site1, site2):
        # Initialize a new sparse NxN matrix
        U = sp.sparse.lil_matrix( (self.N,self.N), dtype=np.complex_ )

        # Loop over all the states
        for state in self.state_to_index.keys():
            from_index = self.state_to_index[state]

            newState = [int(a) for a in state]

            if( state[site1] == "0" and state[site2] == "1" ):
                newState[site1] = 1; newState[site2] = 0;
            elif( state[site1] == "1" and state[site2] == "0" ):
                newState[site1] = 0; newState[site2] = 1;
            else:
                U[from_index, from_index] = 1; # Diagonal untouched
                continue # Skip the rest

            to_index = self.state_to_index[ "".join([str(a) for a in newState])]
            U[to_index, from_index] = 1j

        return U.tocsr()
    
    def constructSqrtSwap(self, site1, site2):
        # Initialize a new sparse NxN matrix
        U = sp.sparse.lil_matrix( (self.N,self.N), dtype=np.complex_ )

        # Loop over all the states
        for state in self.state_to_index.keys():
            from_index = self.state_to_index[state]

            newState = [int(a) for a in state]
            newState[site1], newState[site2] = state[site2], state[site1]
            to_index = self.state_to_index[ "".join([str(a) for a in newState])]
            
            if( to_index == from_index ):
                U[from_index, from_index] = 1; # Diagonal untouched
            else:
                U[from_index,from_index] = 0.5*(1+1j)
                U[to_index, to_index] = 0.5*(1+1j)
                U[to_index, from_index] = 0.5*(1-1j)

        return U.tocsr()

    def applyUnitary( self, U ):
        self.state = U.dot( self.state )
        
    def extend( self, spaces ):
        # add the other space to the front
        # so that all states (0012, 1010 -> 00012, 01010, 10012, 11010)
        
        # Update the localHilbertspaces
        self.localHilbertSizes = list(spaces) + list(self.localHilbertSizes)
        self.N = np.product(self.localHilbertSizes)

        # Build fast indexing dictionariesself.state_to_index = {}
        self.state_to_index = {}
        self.index_to_state = {}
        for i in range(self.size):
            # Convert this integer to the respective bases
            mixedBaseRepresentation = convertToMixedBase(i, self.localHilbertSizes)
            self.state_to_index[mixedBaseRepresentation] = i;
            self.index_to_state[i] = mixedBaseRepresentation;
            
        # Update the state
        newState = sp.sparse.lil_matrix( (self.N,1), dtype=np.complex_ )
        
        for i in self.state.nonzero():
            #r,c = i[0],i[1]
            newState[i] = self.state[i]
            
        self.state = newState

    @property
    def size(self):
        return np.product(self.localHilbertSizes)